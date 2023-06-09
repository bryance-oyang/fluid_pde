/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include <cstdio>

#include "util.hh"
#include "grid.hh"

static inline number vl_lim(number r)
{
	number fabsr;

	fabsr = fabs(r);
	return (r + fabsr) / (1 + fabsr);
}

static inline void plm(number *pql, number *pqr, number q1, number q2, number q3)
{
	number ql, qr, half_step;

	if (q3 - q2 == 0) {
		half_step = 0;
	} else {
		half_step = 0.5 * vl_lim((q2 - q1) / (q3 - q2)) * (q3 - q2);
	}

	ql = q2 - half_step;
	qr = q2 + half_step;

	*pql = ql;
	*pqr = qr;
}

static inline void ppm_lim_parabola(number *pql, number *pqr, number q0, number q1, number q2, number q3, number q4)
{
	number ql = *pql;
	number qr = *pqr;
	number curvl, curvr, curvc, curvf, curv;
	number D;

	D = 1.26;
	if (PPM_STRICT_LIM) {
		D = 1;
	}

	// at local extrema
	if ((qr - q2) * (q2 - ql) <= 0 || (q3 - q2) * (q2 - q1) <= 0) {
		curvc = (q1 + q3) - 2*q2;
		curvl = (q0 + q2) - 2*q1;
		curvr = (q2 + q4) - 2*q3;
		if (WEIRD_PPM) {
			curvf = 4*((ql + qr) - 2*q2);
		} else {
			curvf = 6*((ql + qr) - 2*q2);
		}
		if (SIGN(curvl) == SIGN(curvc) && SIGN(curvc) == SIGN(curvr) && SIGN(curvc) == SIGN(curvf)) {
			curv = SIGN(curvf) * util::fmin4(D*fabs(curvl), D*fabs(curvc), D*fabs(curvr), fabs(curvf));
		} else {
			curv = 0;
		}

		// choose smoothly between q2 and ql/qr
		if (curvf != 0) {
			ql = q2 + (ql - q2) * curv / curvf;
			qr = q2 + (qr - q2) * curv / curvf;
		} else {
			ql = q2;
			qr = q2;
		}
	} else if (fabs(ql - q2) >= 2*fabs(qr - q2)) {
		// move parabola peak out of cell
		ql = q2 - 2*(qr - q2);
	} else if (fabs(qr - q2) >= 2*fabs(ql - q2)) {
		// move parabola peak out of cell
		qr = q2 - 2*(ql - q2);
	}

	*pql = ql;
	*pqr = qr;
}

static inline void fancy_ppm(number *pql, number *pqr, number q0, number q1, number q2, number q3, number q4)
{
	number ql, qr;
	number curvl, curvr, curvf, curv;
	number C;

	C = 1.26;
	if (PPM_STRICT_LIM) {
		C = 1;
	}

	// ql and limiting
	ql = (7.0*(q1 + q2) - (q0 + q3)) / 12;
	curvl = (q0 + q2) - 2*q1;
	curvr = (q1 + q3) - 2*q2;
	curvf = 3*((q1 + q2) - 2*ql);
	// if not monotonic
	if (PPM_ALWAYS_LIM || (curvr - curvf) * (curvl - curvf) > 0) {
		if (SIGN(curvl) == SIGN(curvf) && SIGN(curvf) == SIGN(curvr)) {
			// smooth but potentially apply limiting
			curv = SIGN(curvf) * util::fmin3(C*fabs(curvl), C*fabs(curvr), fabs(curvf));
		} else {
			// 2nd deriv varies a lot
			curv = 0;
		}
		ql = 0.5 * (q1 + q2) - curv / 6;
	}

	// qr and limiting
	qr = (7.0*(q2 + q3) - (q1 + q4)) / 12;
	curvl = (q1 + q3) - 2*q2;
	curvr = (q2 + q4) - 2*q3;
	curvf = 3*((q2 + q3) - 2*qr);
	// if not monotonic
	if (PPM_ALWAYS_LIM || (curvr - curvf) * (curvl - curvf) > 0) {
		if (SIGN(curvl) == SIGN(curvf) && SIGN(curvf) == SIGN(curvr)) {
			// smooth but potentially apply limiting
			curv = SIGN(curvf) * util::fmin3(C*fabs(curvl), C*fabs(curvr), fabs(curvf));
		} else {
			// 2nd deriv varies a lot
			curv = 0;
		}
		qr = 0.5 * (q2 + q3) - curv / 6;
	}

	*pql = ql;
	*pqr = qr;

	ppm_lim_parabola(pql, pqr, q0, q1, q2, q3, q4);

	if (PPM_STRICT_LIM) {
		*pql = fmin(fmax(q1,q2),*pql);
		*pql = fmax(fmin(q1,q2),*pql);

		*pqr = fmin(fmax(q3,q2),*pqr);
		*pqr = fmax(fmin(q3,q2),*pqr);
	}
}

void Grid::Reconstruct(int dir)
{
	int di, dj;
	if (dir == 0) {
		di = 1;
		dj = 0;
	} else {
		di = 0;
		dj = 1;
	}

	for (int m = 0; m < NQUANT; m++) {
		for (int i = ilr; i < iur; i++) {
			// cell loop
			for (int j = jl-1; j < ju+1; j++) {
				number ql, qr;
				number q0 = prim(m,i-2*di,j-2*dj);
				number q1 = prim(m,i-di,j-dj);
				number q2 = prim(m,i,j);
				number q3 = prim(m,i+di,j+dj);
				number q4 = prim(m,i+2*di,j+2*dj);

				if (reconstruct_order == 1) {
					ql = q2;
					qr = q2;
				} else if (reconstruct_order == 2) {
					plm(&ql, &qr, q1, q2, q3);
				} else if (reconstruct_order == 3) {
					fancy_ppm(&ql, &qr, q0, q1, q2, q3, q4);
				} else {
					ql = 0;
					qr = 0;
				}

				Rprim(m,i,j) = ql;
				Lprim(m,i+di,j+dj) = qr;
			}
		}
	}
}
