/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "riemann.hh"
#include "macro.hh"

namespace riemann {

void HLLC(const Array<number> &Lprim, const Array<number> &Lcons, const Array<number> &Lw_array,
	const Array<number> &Rprim, const Array<number> &Rcons, const Array<number> &Rw_array,
	Array<number> &J, int dir, int il, int iuf, int jl, int ju)
{
	for (int i = il; i < iuf; i++) {
		//face loop
		for (int j = jl; j < ju+1; j++) {
			int m;
			number Lw, Rw, Mw;
			number Lv, Rv;
			number Lpress, Rpress, Mpress;
			number Lrho, Rrho, rho2, rho3;
			number Le, Re, e2, e3;

			Lw = Lw_array(i,j);
			Rw = Rw_array(i,j);

			if (Lw == 0 && Rw == 0) {
				for (m = 0; m < NQUANT; m++) {
					J(m,i,j) = 0;
				}
				continue;
			}

			Lrho = Lprim(0,i,j);
			Rrho = Rprim(0,i,j);
			Lv = Lprim(1+dir,i,j);
			Rv = Rprim(1+dir,i,j);
			Lpress = Lprim(3,i,j);
			Rpress = Rprim(3,i,j);
			Le = Lcons(3,i,j);
			Re = Rcons(3,i,j);

			// supersonic
			if (Rw < 0) {
				for (m = 0; m < NQUANT; m++) {
					J(m,i,j) = Rcons(m,i,j) * Rv;
					if (m == 1+dir) {
						J(m,i,j) += Rpress;
					}
					if (m == 3) {
						J(m,i,j) += Rpress * Rv;
					}
				}
				continue;
			}
			if (Lw > 0) {
				for (m = 0; m < NQUANT; m++) {
					J(m,i,j) = Lcons(m,i,j) * Lv;
					if (m == 1+dir) {
						J(m,i,j) += Lpress;
					}
					if (m == 3) {
						J(m,i,j) += Lpress * Lv;
					}
				}
				continue;
			}

			// middle wave and intermediate left/right density
			Mw = ((Rrho*Rv*(Rv-Rw) + Rpress) - (Lrho*Lv*(Lv-Lw) + Lpress)) / (Rrho*(Rv-Rw) - Lrho*(Lv-Lw));
			rho2 = Lrho * (Lv - Lw) / (Mw - Lw);
			rho3 = Rrho * (Rv - Rw) / (Mw - Rw);

			if (Mw > 0) {
				Mpress = Lrho*SQR(Lv) + Lpress - Lw*Lrho*Lv - rho2*SQR(Mw) + Lw*rho2*Mw;
			} else if (Mw < 0) {
				Mpress = Rrho*SQR(Rv) + Rpress - Rw*Rrho*Rv - rho3*SQR(Mw) + Rw*rho3*Mw;
			} else {
				Mpress = 0.5 * ((Lrho*SQR(Lv) + Lpress - Lw*Lrho*Lv - rho2*SQR(Mw) + Lw*rho2*Mw)
						+ (Rrho*SQR(Rv) + Rpress - Rw*Rrho*Rv - rho3*SQR(Mw) + Rw*rho3*Mw));
			}

			// contact wave in middle
			if (Mw == 0) {
				for (m = 0; m < NQUANT; m++) {
					J(m,i,j) = 0;
					if (m == 1+dir) {
						J(m,i,j) += Mpress;
					}
				}
				continue;
			}

			if (Mw < 0) {
				// intermediate energy
				e3 = (Rv*(Re+Rpress) - Rw*Re - Mw*Mpress) / (Mw - Rw);

				J(0,i,j) = rho3 * Mw;
				if (dir == 0) {
					J(1,i,j) = rho3 * SQR(Mw) + Mpress;
					J(2,i,j) = rho3 * Rprim(2,i,j) * Mw;
				} else {
					J(1,i,j) = rho3 * Rprim(1,i,j) * Mw;
					J(2,i,j) = rho3 * SQR(Mw) + Mpress;
				}
				J(3,i,j) = (e3 + Mpress) * Mw;
				for (m = 4; m < NQUANT; m++) {
					J(m,i,j) = rho3 * Rprim(m,i,j) * Mw;
				}
			} else {
				// intermediate energy
				e2 = (Lv*(Le+Lpress) - Lw*Le - Mw*Mpress) / (Mw - Lw);

				J(0,i,j) = rho2 * Mw;
				if (dir == 0) {
					J(1,i,j) = rho2 * SQR(Mw) + Mpress;
					J(2,i,j) = rho2 * Lprim(2,i,j) * Mw;
				} else {
					J(1,i,j) = rho2 * Lprim(1,i,j) * Mw;
					J(2,i,j) = rho2 * SQR(Mw) + Mpress;
				}
				J(3,i,j) = (e2 + Mpress) * Mw;
				for (m = 4; m < NQUANT; m++) {
					J(m,i,j) = rho2 * Lprim(m,i,j) * Mw;
				}
			}
		}
	}
}

void HLLE(const Array<number> &Lcons, const Array<number> &LJ_array, const Array<number> &Lw_array,
	const Array<number> &Rcons, const Array<number> &RJ_array, const Array<number> &Rw_array,
	Array<number> &J, int dir, int il, int iuf, int jl, int ju)
{
	(void)dir;

	for (int m = 0; m < NQUANT; m++) {
		for (int i = il; i < iuf; i++) {
			//face loop
			for (int j = jl; j < ju+1; j++) {
				number Lq, Rq, LJ, RJ, Lw, Rw;

				Lq = Lcons(m,i,j);
				Rq = Rcons(m,i,j);
				LJ = LJ_array(m,i,j);
				RJ = RJ_array(m,i,j);
				Lw = Lw_array(i,j);
				Rw = Rw_array(i,j);

				if (Lw == 0 && Rw == 0) {
					J(m,i,j) = 0;
				} else if (Rw <= 0) {
					J(m,i,j) = RJ;
				} else if (Lw >= 0) {
					J(m,i,j) = LJ;
				} else {
					J(m,i,j) = (LJ*Rw - RJ*Lw + Rw*Lw*(Rq - Lq)) / (Rw - Lw);
				}
			}
		}
	}
}

} // namespace riemann
