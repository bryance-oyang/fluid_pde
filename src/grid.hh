/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef GRID_H
#define GRID_H

#include <mutex>
#include <condition_variable>
#include "array.hh"
#include "macro.hh"
#include "config.hh"

#define NQUANT ((int)(4+(int)(NSCALAR)))

class Grid {
public:
	int nu;
	int nv;
	number umin;
	number umax;
	number vmin;
	number vmax;
	number du;
	number dv;

	int tid;
	// for nonghost of thread division
	int il;
	int iu;
	int iuf; // for face loops
	int ilr; // for reconstruction
	int iur; // for reconstruction

	int jl;
	int ju;

	int reconstruct_order;

	number &time;
	number &dt;
	// time in middle of rk integration
	number &step_time;
	// dt in middle of rk integration
	number &step_dt;

	number rho_floor;
	number press_floor;
	number gamma;

	Array<number> u_cc;
	Array<number> v_cc;
	Array<number> u_ufc;
	Array<number> v_ufc;
	Array<number> u_vfc;
	Array<number> v_vfc;

	Array<number> cons;
	Array<number> prim;
	Array<number> cons_gen;

	Array<number> fluxdiv;
	Array<number> src;

	Array<number> Ju;
	Array<number> Jv;

	// reconstruction vars
	Array<number> Lprim;
	Array<number> Lcons;
	Array<number> Rprim;
	Array<number> Rcons;

	// wavespeed
	Array<number> Lw;
	Array<number> Rw;

	Grid(number &time, number &dt, number &step_time, number &step_dt);

	// set grid properties
	void Property();

	// called after Property
	void InitGrid();

	// setup
	void AllocGrid();
	void AttachReference(Grid &g);
	void DetachReference();
	void InitUVCoord();

	// setup prim
	void InitCond();

	void ConsLim();
	void ConsToPrim();
	void PointPrimToCons(const Array<number> &prim, Array<number> &cons);
	void PrimLim(Array<number> &prim);
	void PrimToCons(const Array<number> &prim, Array<number> &cons);

	void Reconstruct(int dir);
	void Wavespeed(int dir);
	void CalculateSrc();
	void DetermineDt(int dir);
	void CalculateFluxDiv();

	// boundary
	void Boundary(number time);

	void InflowBoundaryLeft(number rho, number vx, number vy, number press);
	void InflowBoundaryRight(number rho, number vx, number vy, number press);
	void InflowBoundaryBot(number rho, number vx, number vy, number press);
	void InflowBoundaryTop(number rho, number vx, number vy, number press);

	void PeriodicBoundaryLeft();
	void PeriodicBoundaryRight();
	void PeriodicBoundaryBot();
	void PeriodicBoundaryTop();

	void PeriodicBoundaryLB();
	void PeriodicBoundaryRB();
	void PeriodicBoundaryRT();
	void PeriodicBoundaryLT();

	void SmoothBoundaryLeft();
	void SmoothBoundaryRight();
	void SmoothBoundaryBot();
	void SmoothBoundaryTop();

	void SmoothBoundaryLB();
	void SmoothBoundaryRB();
	void SmoothBoundaryRT();
	void SmoothBoundaryLT();

	void ReflectingBoundaryLeft();
	void ReflectingBoundaryRight();
	void ReflectingBoundaryBot();
	void ReflectingBoundaryTop();

	void ReflectingBoundaryLB();
	void ReflectingBoundaryRB();
	void ReflectingBoundaryRT();
	void ReflectingBoundaryLT();
};

#endif /* GRID_H */
