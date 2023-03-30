/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef INIT_COND_H
#define INIT_COND_H

#include "../config.hh"
#include "../integrator.hh"
#include "../grid.hh"

#include <math.h>

void Integrator::Property()
{
	// max steps
	max_epoch = 1000000000;
	// max output
	max_out = 50;
	// max output time
	out_tf = 1;

	// set cfl number
	cfl_num = 0.43;

	// time integrator choice
	//Euler();
	//RK2();
	//SSPRK3();
	SSPRK4();

	// autocompute
	out_dt = out_tf / (max_out - 1);
}

void Grid::Property()
{
	// 1: 1st order no reconstruction, 2: 2nd order linear, 3: 4th order parabolic
	reconstruct_order = 3;

	// floors
	rho_floor = 1e-8;
	press_floor = 1e-10;

	// grid resolution
	nu = NU;
	nv = NV;

	// coordinate bound
	umin = -1;
	umax = 1;
	vmin = -1;
	vmax = 1;

	// adiabatic index
	gamma = 1.4;
}

// initial conditions of grid
void Grid::InitCond()
{
	Array<number> tmp_prim(NQUANT);
	Array<number> tmp_cons(NQUANT);

	for (int i = 0; i < nu; i++) {
		for (int j = 0; j < nv; j++) {
			// coordinates (x, y)
			number x, y;
			x = u_cc(i);
			y = v_cc(j);

			// 0: density; 1: vel1; 2: vel2; 3: pressure, 4-: passive scalars (set num in config.h)
			tmp_prim(0) = 1;
			tmp_prim(1) = 0;
			tmp_prim(2) = 0;
			tmp_prim(3) = 1;

			PointPrimToCons(tmp_prim, tmp_cons);

			cons_gen(0,i,j) = tmp_cons(0);
			cons_gen(1,i,j) = tmp_cons(1);
			cons_gen(2,i,j) = tmp_cons(2);
			cons_gen(3,i,j) = tmp_cons(3);
		}
	}

	cons.copy_data_from(cons_gen);
}

void Grid::CalculateSrc()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = il; i < iu; i++) {
			for (int j = jl; j < ju; j++) {
				src(m,i,j) = 0;
			}
		}
	}
}

void Grid::Boundary(number time)
{
	//edges
	PeriodicBoundaryLeft();
	PeriodicBoundaryRight();
	PeriodicBoundaryBot();
	PeriodicBoundaryTop();
	// corners
	PeriodicBoundaryLB();
	PeriodicBoundaryRB();
	PeriodicBoundaryRT();
	PeriodicBoundaryLT();
}

#endif /* INIT_COND_H */
