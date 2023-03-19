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

#include "ws_ctube.h"
#include "grid.h"

class GridConverter {
public:
	Array<double> preimage;
	Array<uint8_t> image;

	GridConverter() {
		preimage = Array<double>{NU, NV};
		image = Array<uint8_t>{NU, NV, 3};
	}

	double clip(double x)
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

class Broadcast {
public:
	std::unique_ptr<std::thread> thread;
	ws_ctube *ctube = NULL;
	std::atomic<int> should_terminate;

	Grid &g;
	GridConverter converter;

	Broadcast(Grid &g, int port, int max_nclient, int timeout_ms, double max_broadcast_fps)
	: g{g} {
		if (start_ctube(port, max_nclient, timeout_ms, max_broadcast_fps)) {
			start_thread();
		}
	}

	~Broadcast() noexcept {
		join();
	}

	bool start_ctube(int port, int max_nclient, int timeout_ms, double max_broadcast_fps)
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

	bool start_thread()
	{
		stop_thread();
		should_terminate.store(0);
		thread = std::make_unique<std::thread>(&Broadcast::thread_main, this);
		return true;
	}
	void stop_thread()
	{
		if (thread) {
			should_terminate.store(1);
			if (thread->joinable()) {
				thread->join();
			}
			thread.reset();
		}
	}

	void join()
	{
		stop_thread();
		stop_ctube();
	}

	void thread_main()
	{
		for (;;) {
			{ /* lock camera mutex */
				std::unique_lock<std::mutex> mutex{g.mutex};
				while (!g.broadcast_signal) {
					using namespace std::chrono_literals;
					g.cond.wait_for(mutex, 200ms);

					/* check if we should exit */
					if (should_terminate.load()) {
						return;
					}
				}
				g.broadcast_signal = false;
				converter.make_image(g);
			} /* unlock camera mutex */

			ws_ctube_broadcast(ctube, converter.image.data,
				converter.image.bytes());
		}
	}
};

#endif /* BROADCAST_H */
