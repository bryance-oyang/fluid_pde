/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "grid.h"
#include "riemann.h"

#include <math.h>
#include <cstdio>
#include <stdlib.h>

Grid::Grid(double &time, double &dt, double &step_time, double &step_dt)
: time{time}, dt{dt}, step_time{step_time}, step_dt{step_dt}
{
	Property();
}

void Grid::InitGrid()
{
	if (reconstruct_order != 1 && reconstruct_order != 2 && reconstruct_order != 3) {
		printf("Bad reconstruct_order\n");
		exit(EXIT_FAILURE);
	}

	AllocGrid();
	InitUVCoord();

	InitCond();

	ConsLim();
	ConsToPrim();
	Boundary(time);
	ConsLim();
	ConsToPrim();
}

void Grid::AllocGrid()
{
	// include boundary
	nu += 2*NGHOST;
	nv += 2*NGHOST;

	// calculated without boundary
	du = (umax - umin) / (nu - 2*NGHOST);
	dv = (vmax - vmin) / (nv - 2*NGHOST);

	time = 0;
	dt = 0;

	// coordinates
	u_cc = Array<double>{nu};
	v_cc = Array<double>{nv};
	u_ufc = Array<double>{nu+1};
	v_ufc = Array<double>{nv+1};
	u_vfc = Array<double>{nu+1};
	v_vfc = Array<double>{nv+1};

	// hydro
	cons = Array<double>{NQUANT, nu, nv};
	prim = Array<double>{NQUANT, nu, nv};
	cons_gen = Array<double>{NQUANT, nu, nv};

	fluxdiv = Array<double>{NQUANT, nu, nv};
	src = Array<double>{NQUANT, nu, nv};

	// current
	Ju = Array<double>{NQUANT, nu+1, nv+1};
	Jv = Array<double>{NQUANT, nu+1, nv+1};

	// reconstruction vars
	Lprim = Array<double>{NQUANT, nu+1, nv+1};
	Lcons = Array<double>{NQUANT, nu+1, nv+1};
	Rprim = Array<double>{NQUANT, nu+1, nv+1};
	Rcons = Array<double>{NQUANT, nu+1, nv+1};

	// wavespeed
	Lw = Array<double>{nu+1, nv+1};
	Rw = Array<double>{nu+1, nv+1};
}

void Grid::AttachReference(Grid &g)
{
	nu = g.nu;
	nv = g.nv;
	umin = g.umin;
	umax = g.umax;
	vmin = g.vmin;
	vmax = g.vmax;
	du = g.du;
	dv = g.dv;

	u_cc.attach_reference(g.u_cc);
	v_cc.attach_reference(g.v_cc);
	u_ufc.attach_reference(g.u_ufc);
	v_ufc.attach_reference(g.v_ufc);
	u_vfc.attach_reference(g.u_vfc);
	v_vfc.attach_reference(g.v_vfc);
	cons.attach_reference(g.cons);
	prim.attach_reference(g.prim);
	cons_gen.attach_reference(g.cons_gen);
	fluxdiv.attach_reference(g.fluxdiv);
	src.attach_reference(g.src);
	Ju.attach_reference(g.Ju);
	Jv.attach_reference(g.Jv);
	Lprim.attach_reference(g.Lprim);
	Lcons.attach_reference(g.Lcons);
	Rprim.attach_reference(g.Rprim);
	Rcons.attach_reference(g.Rcons);
	Lw.attach_reference(g.Lw);
	Rw.attach_reference(g.Rw);

	reconstruct_order = g.reconstruct_order;
	rho_floor = g.rho_floor;
	press_floor = g.press_floor;
	gamma = g.gamma;
}

void Grid::DetachReference()
{
	u_cc.detach_reference();
	v_cc.detach_reference();
	u_ufc.detach_reference();
	v_ufc.detach_reference();
	u_vfc.detach_reference();
	v_vfc.detach_reference();
	cons.detach_reference();
	prim.detach_reference();
	cons_gen.detach_reference();
	fluxdiv.detach_reference();
	src.detach_reference();
	Ju.detach_reference();
	Jv.detach_reference();
	Lprim.detach_reference();
	Lcons.detach_reference();
	Rprim.detach_reference();
	Rcons.detach_reference();
	Lw.detach_reference();
	Rw.detach_reference();
}

void Grid::InitUVCoord()
{
	for (int i = 0; i < nu; i++) {
		for (int j = 0; j < nv; j++) {
			u_cc(i) = umin + (i-NGHOST+0.5)*du;
			v_cc(j) = vmin + (j-NGHOST+0.5)*dv;
		}
	}

	for (int i = 0; i < nu+1; i++) {
		for (int j = 0; j < nv+1; j++) {
			u_ufc(i) = umin + (i-NGHOST+0.5)*du;
			v_ufc(j) = vmin + (j-NGHOST)*dv;

			u_vfc(i) = umin + (i-NGHOST)*du;
			v_vfc(j) = vmin + (j-NGHOST+0.5)*dv;
		}
	}
}
