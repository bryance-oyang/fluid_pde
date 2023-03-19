/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <math.h>

#include "grid.h"
#include "riemann.h"

void Grid::ConsLim()
{
	ConsToPrim();
	PrimLim(prim);
	PrimToCons(prim, cons);
}

void Grid::PrimLim(Array<double> &prim)
{
	for (int i = 0; i < prim.n[1]; i++) {
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

void Grid::PrimToCons(const Array<double> &prim, Array<double> &cons)
{
	for (int i = 0; i < prim.n[1]; i++) {
		for (int j = 0; j < prim.n[2]; j++) {
			double rho = prim(0,i,j);
			double vsquared = SQR(prim(1,i,j)) + SQR(prim(2,i,j));

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
	for (int i = 0; i < cons.n[1]; i++) {
		for (int j = 0; j < cons.n[2]; j++) {
			double rho = cons(0,i,j);
			double v1 = cons(1,i,j) / rho;
			double v2 = cons(2,i,j) / rho;
			double ke = 0.5 * rho * (SQR(v1) + SQR(v2));

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

void Grid::PointPrimToCons(const Array<double> &prim, Array<double> &cons)
{
	double rho = prim(0);
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
	int di, dj, pdi, pdj;
	if (dir == 0) {
		di = 1;
		dj = 0;
		pdi = 0;
		pdj = 1;
	} else {
		di = 0;
		dj = 1;
		pdi = 1;
		pdj = 0;
	}

	for (int i = NGHOST-pdi; i < nu-NGHOST+1; i++) {
		// face loop
		for (int j = NGHOST-pdj; j < nv-NGHOST+1; j++) {
			double Lcs, Rcs, Lv, Rv;

			Lcs = fmax(sqrt(gamma * Lprim(3,i,j) / Lprim(0,i,j)), sqrt(gamma * prim(3,i-di,j-dj) / prim(0,i-di,j-dj)));
			Rcs = fmax(sqrt(gamma * Rprim(3,i,j) / Rprim(0,i,j)), sqrt(gamma * prim(3,i,j) / prim(0,i,j)));
			Lv = fmin(Lprim(1+dir,i,j), prim(1+dir,i-di,j-dj));
			Rv = fmax(Rprim(1+dir,i,j), prim(1+dir,i,j));

			Lw(i,j) = fmin(Lv - Lcs, Rv - Rcs);
			Rw(i,j) = fmax(Lv + Lcs, Rv + Rcs);
		}
	}
}

void Grid::CalculateRiemannJ(int dir)
{
	Array<double> *J;
	if (dir == 0) {
		J = &Ju;
	} else {
		J = &Jv;
	}

	Reconstruct(dir);

	PrimLim(Lprim);
	PrimLim(Rprim);

	PrimToCons(Lprim, Lcons);
	PrimToCons(Rprim, Rcons);

	Wavespeed(dir);

	riemann::HLLC(Lprim, Lcons, Lw, Rprim, Rcons, Rw, *J, dir);
}

void Grid::CalculateJ(Array<double> &J, int dir)
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < J.n[1]; i++) {
			for (int j = 0; j < J.n[2]; j++) {
				J(m,i,j) = prim(1+dir,i,j) * cons(m,i,j);

				if (m == 1+dir) {
					J(m,i,j) += prim(3,i,j);
				}
				if (m == 3) {
					J(m,i,j) += prim(1+dir,i,j) * prim(3,i,j);
				}
			}
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
	double ds;
	if (dir == 0) {
		di = 1;
		dj = 0;
		ds = du;
	} else {
		di = 0;
		dj = 1;
		ds = dv;
	}

	for (int i = NGHOST; i < nu-NGHOST; i++) {
		// cell loop
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			double w1, w2;
			double cross_time;

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
		for (int i = NGHOST; i < nu-NGHOST; i++) {
			// cell loop
			for (int j = NGHOST; j < nv-NGHOST; j++) {
					fluxdiv(m,i,j) = (Ju(m,i,j) - Ju(m,i+1,j)) / du
						+ (Jv(m,i,j) - Jv(m,i,j+1)) / dv;
			}
		}
	}
}
