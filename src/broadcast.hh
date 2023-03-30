/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef BROADCAST_H
#define BROADCAST_H

#include <atomic>
#include <thread>
#include <mutex>
#include <chrono>
#include <cstdint>
#include <cmath>

#include "ws_ctube.hh"
#include "grid.hh"

class GridConverter {
public:
	Array<number> preimage;
	Array<uint8_t> image;

	GridConverter() {
		preimage = Array<number>{NU, NV};
		image = Array<uint8_t>{NU, NV, 3};
	}

	number clip(number x)
	{
		return fmax(BROADCAST_PREIMAGE_MIN, fmin(BROADCAST_PREIMAGE_MAX, x));
	}

	void make_image(const Grid &g)
	{
		for (int i = 0; i < preimage.n[0]; i++) {
			for (int j = 0; j < preimage.n[1]; j++) {
				preimage(i, j) = clip(log10(g.cons(0, i+NGHOST, j+NGHOST)));
			}
		}
		image_from_preimage();
	}

	void image_from_preimage()
	{
		for (int i = 0; i < preimage.n[0]; i++) {
			for (int j = 0; j < preimage.n[1]; j++) {
				uint8_t value = (uint8_t)(255.001 * (preimage(i, j) - BROADCAST_PREIMAGE_MIN) / (BROADCAST_PREIMAGE_MAX - BROADCAST_PREIMAGE_MIN));
				for (int k = 0; k < 3; k++) {
					image(i, j, k) = value;
				}
			}
		}
	}
};

class Broadcaster {
public:
	ws_ctube *ctube = NULL;
	Grid &g;
	GridConverter converter;

	Broadcaster(Grid &g, int port, int max_nclient, int timeout_ms, number max_broadcast_fps)
	: g{g} {
		start_ctube(port, max_nclient, timeout_ms, max_broadcast_fps);
	}

	bool start_ctube(int port, int max_nclient, int timeout_ms, number max_broadcast_fps)
	{
		stop_ctube();
		ctube = ws_ctube_open(port, max_nclient, timeout_ms, max_broadcast_fps);
		return ctube != NULL;
	}
	void stop_ctube()
	{
		if (ctube != NULL) {
			ws_ctube_close(ctube);
			ctube = NULL;
		}
	}

	void broadcast() {
		converter.make_image(g);
		ws_ctube_broadcast(ctube, converter.image.data,
			converter.image.bytes());
	}
};

#endif /* BROADCAST_H */
