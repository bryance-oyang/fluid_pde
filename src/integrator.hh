/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "grid.hh"

class Integrator {
public:
	int s;
	int nstep;
	Array<double> weight;
	Array<double> time_weight;

	double cfl_num;

	int max_epoch;
	int max_out;
	double out_tf;
	double out_dt;

	// ssprk4
	/*
	int ssprk4;
	Array rk4_fin_weight;
	Array rk4_u2;
	Array rk4_u3;
	Array rk4_deriv3;
	*/

	void Property();

	void ComputeTimeWeight();

	void Euler();
	void RK2();
	void SSPRK3();
	//void SSPRK4();

	void AddFluxDivSrc(Grid *g);
};

#endif /* INTEGRATOR_H */
