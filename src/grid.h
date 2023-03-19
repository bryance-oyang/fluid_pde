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
#include "array.h"
#include "macro.h"
#include "config.h"

#define NQUANT ((int)(4+(int)(NSCALAR)))

class Grid {
public:
	std::mutex mutex;
	std::condition_variable cond;
	bool broadcast_signal = false;

	int nu;
	int nv;
	double umin;
	double umax;
	double vmin;
	double vmax;
	double du;
	double dv;

	int il;
	int iu;
	int jl;
	int ju;

	int reconstruct_order;

	double &time;
	double &dt;
	// time in middle of rk integration
	double &step_time;
	// dt in middle of rk integration
	double &step_dt;

	double rho_floor;
	double press_floor;
	double gamma;

	Array<double> u_cc;
	Array<double> v_cc;
	Array<double> u_ufc;
	Array<double> v_ufc;
	Array<double> u_vfc;
	Array<double> v_vfc;

	Array<double> cons;
	Array<double> prim;
	Array<double> cons_gen;

	Array<double> fluxdiv;
	Array<double> src;

	Array<double> Ju;
	Array<double> Jv;

	// reconstruction vars
	Array<double> Lprim;
	Array<double> Lcons;
	Array<double> Rprim;
	Array<double> Rcons;

	// wavespeed
	Array<double> Lw;
	Array<double> Rw;

	Grid(double &time, double &dt, double &step_time, double &step_dt);

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
	void PointPrimToCons(const Array<double> &prim, Array<double> &cons);
	void PrimLim(Array<double> &prim);
	void PrimToCons(const Array<double> &prim, Array<double> &cons);

	void Reconstruct(int dir);
	void Wavespeed(int dir);
	void CalculateJ(Array<double> &J, int dir);
	void CalculateSrc();
	void DetermineDt(int dir);
	void CalculateFluxDiv();

	// boundary
	void Boundary(double time);

	void InflowBoundaryLeft(double rho, double vx, double vy, double press);
	void InflowBoundaryRight(double rho, double vx, double vy, double press);
	void InflowBoundaryBot(double rho, double vx, double vy, double press);
	void InflowBoundaryTop(double rho, double vx, double vy, double press);

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
