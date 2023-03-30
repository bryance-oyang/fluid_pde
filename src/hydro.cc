/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cmath>

#include "grid.hh"
#include "riemann.hh"

// for [0, lim)
static void determine_loop_limits(int tid, int lim, int *iil, int *iiu)
{
	if (tid >= 0) {
		int nii = (lim + NTHREAD - 1) / NTHREAD;
		*iil = tid * nii;
		*iiu = std::min((tid + 1)*nii, lim);
	} else {
		// global grid
		*iil = 0;
		*iiu = lim;
	}
}

void Grid::ConsLim()
{
	ConsToPrim();
	PrimLim(prim);
	PrimToCons(prim, cons);
}

void Grid::PrimLim(Array<number> &prim)
{
	int iil, iiu;

	determine_loop_limits(tid, prim.n[1], &iil, &iiu);
	for (int i = iil; i < iiu; i++) {
		for (int j = 0; j < prim.n[2]; j++) {
			if (prim(0,i,j) < rho_floor) {
				prim(0,i,j) = rho_floor;
			}

			if (prim(3,i,j) < press_floor) {
				prim(3,i,j) = press_floor;
			}

			for (int m = 4; m < NQUANT; m++) {
				if (prim(m,i,j) < 0) {
					prim(m,i,j) = 0;
				}
			}
		}
	}
}

void Grid::PrimToCons(const Array<number> &prim, Array<number> &cons)
{
	int iil, iiu;

	determine_loop_limits(tid, prim.n[1], &iil, &iiu);
	for (int i = iil; i < iiu; i++) {
		for (int j = 0; j < prim.n[2]; j++) {
			number rho = prim(0,i,j);
			number vsquared = SQR(prim(1,i,j)) + SQR(prim(2,i,j));

			cons(0,i,j) = rho;
			cons(1,i,j) = rho * prim(1,i,j);
			cons(2,i,j) = rho * prim(2,i,j);
			cons(3,i,j) = 0.5 * rho * vsquared + prim(3,i,j) / (gamma - 1);

			for (int m = 4; m < NQUANT; m++) {
				cons(m,i,j) = rho * prim(m,i,j);
			}
		}
	}
}

void Grid::ConsToPrim()
{
	int iil, iiu;

	determine_loop_limits(tid, cons.n[1], &iil, &iiu);
	for (int i = iil; i < iiu; i++) {
		for (int j = 0; j < cons.n[2]; j++) {
			number rho = cons(0,i,j);
			number v1 = cons(1,i,j) / rho;
			number v2 = cons(2,i,j) / rho;
			number ke = 0.5 * rho * (SQR(v1) + SQR(v2));

			prim(0,i,j) = rho;
			prim(1,i,j) = v1;
			prim(2,i,j) = v2;
			prim(3,i,j) = (cons(3,i,j) - ke) * (gamma - 1);

			for (int m = 4; m < NQUANT; m++) {
				prim(m,i,j) = cons(m,i,j) / rho;
			}
		}
	}
}

void Grid::PointPrimToCons(const Array<number> &prim, Array<number> &cons)
{
	number rho = prim(0);
		cons(0) = rho;
		cons(1) = rho * prim(1);
		cons(2) = rho * prim(2);
		cons(3) = 0.5 * rho * (SQR(prim(1)) + SQR(prim(2))) + prim(3) / (gamma - 1);
	for (int m = 4; m < NQUANT; m++) {
		cons(m) = rho * prim(m);
	}
}

void Grid::Wavespeed(int dir)
{
	int di, dj;
	if (dir == 0) {
		di = 1;
		dj = 0;
	} else {
		di = 0;
		dj = 1;
	}

	for (int i = il; i < iuf; i++) {
		// face loop
		for (int j = jl; j < ju+1; j++) {
			number Lcs, Rcs, Lv, Rv;

			Lcs = fmax(sqrt(gamma * Lprim(3,i,j) / Lprim(0,i,j)), sqrt(gamma * prim(3,i-di,j-dj) / prim(0,i-di,j-dj)));
			Rcs = fmax(sqrt(gamma * Rprim(3,i,j) / Rprim(0,i,j)), sqrt(gamma * prim(3,i,j) / prim(0,i,j)));
			Lv = fmin(Lprim(1+dir,i,j), prim(1+dir,i-di,j-dj));
			Rv = fmax(Rprim(1+dir,i,j), prim(1+dir,i,j));

			Lw(i,j) = fmin(Lv - Lcs, Rv - Rcs);
			Rw(i,j) = fmax(Lv + Lcs, Rv + Rcs);
		}
	}
}

void __attribute__((weak)) Grid::CalculateSrc()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < nu; i++) {
			for (int j = 0; j < nv; j++) {
				src(m,i,j) = 0*prim(0)*cons(0)*time*dt;
			}
		}
	}
}

void Grid::DetermineDt(int dir)
{
	int di, dj;
	number ds;
	if (dir == 0) {
		di = 1;
		dj = 0;
		ds = du;
	} else {
		di = 0;
		dj = 1;
		ds = dv;
	}

	// not thread local
	for (int i = NGHOST; i < nu-NGHOST; i++) {
		// cell loop
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			number w1, w2;
			number cross_time;

			w1 = fabs(Rw(i,j));
			w2 = fabs(Lw(i+di,j+dj));

			cross_time = fmin(ds/w1, ds/w2);

			if (cross_time < dt) {
				dt = cross_time;
			}
		}
	}
}

void Grid::CalculateFluxDiv()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = il; i < iu; i++) {
			// cell loop
			for (int j = jl; j < ju; j++) {
				fluxdiv(m,i,j) = (Ju(m,i,j) - Ju(m,i+1,j)) / du
					+ (Jv(m,i,j) - Jv(m,i,j+1)) / dv;
			}
		}
	}
}
