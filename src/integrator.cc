/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "integrator.hh"

#include <cfloat>

void Integrator::ComputeTimeWeight()
{
	time_weight = Array<double>{nstep};

	for (int i = 0; i < nstep; i++) {
		time_weight(i) = (weight(i,1) + weight(i,2)) / (weight(i,0) + weight(i,1));
	}
}

void Integrator::AddFluxDivSrc(Grid *g)
{
	int il = g->il;
	int iu = g->iu;
	int jl = g->jl;
	int ju = g->ju;

	for (int m = 0; m < NQUANT; m++) {
		for (int i = il; i < iu; i++) {
			// cell loop
			for (int j = jl; j < ju; j++) {
				double deriv;

				deriv = g->fluxdiv(m,i,j) + g->src(m,i,j);

				g->cons(m,i,j) = weight(s,0)*g->cons_gen(m,i,j) + weight(s,1)*g->cons(m,i,j) + weight(s,2)*deriv*g->dt;

				// ssprk4 logic
				/*
				if (ssprk4) {
					if (s == 1) {
						rk4_u2(m,i,j) = g->cons(m,i,j);
					}
					if (s == 2) {
						rk4_u3(m,i,j) = g->cons(m,i,j);
					}
					if (s == 3) {
						rk4_deriv3(m,i,j) = deriv;
					}
					if (s == 4) {
						g->cons(m,i,j) += rk4_fin_weight(0) * rk4_u2(m,i,j)
							+ rk4_fin_weight(1) * rk4_u3(m,i,j)
							+ rk4_fin_weight(2) * rk4_deriv3(m,i,j) * g->dt;
					}
				}
				*/
			}
		}
	}
}

void Integrator::Euler()
{
	nstep = 1;
	weight = Array<double>{nstep, 3};

	weight(0,0) = 1;
	weight(0,1) = 0;
	weight(0,2) = 1;

	ComputeTimeWeight();
}

void Integrator::RK2()
{
	nstep = 2;
	weight = Array<double>{nstep, 3};

	weight(0,0) = 1;
	weight(0,1) = 0;
	weight(0,2) = 1;

	weight(1,0) = 1.0/2.0;
	weight(1,1) = 1.0/2.0;
	weight(1,2) = 1.0/2.0;

	ComputeTimeWeight();
}

void Integrator::SSPRK3()
{
	nstep = 3;
	weight = Array<double>{nstep, 3};

	weight(0,0) = 1;
	weight(0,1) = 0;
	weight(0,2) = 1;

	weight(1,0) = 3.0/4.0;
	weight(1,1) = 1.0/4.0;
	weight(1,2) = 1.0/4.0;

	weight(2,0) = 1.0/3.0;
	weight(2,1) = 2.0/3.0;
	weight(2,2) = 2.0/3.0;

	ComputeTimeWeight();
}

/*
void Integrator::SSPRK4()
{
	ssprk4 = 1;
	nstep = 5;
	weight = Array<double>{nstep, 3};

	weight(0,0) = 1;
	weight(0,1) = 0;
	weight(0,2) = 0.391752226571890;

	weight(1,0) = 0.444370493651235;
	weight(1,1) = 0.555629506348765;
	weight(1,2) = 0.368410593050371;

	weight(2,0) = 0.620101851488403;
	weight(2,1) = 0.379898148511597;
	weight(2,2) = 0.251891774271694;

	weight(3,0) = 0.178079954393132;
	weight(3,1) = 0.821920045606868;
	weight(3,2) = 0.544974750228521;

	weight(4,0) = 0;
	weight(4,1) = 0.386708617503269;
	weight(4,2) = 0.226007483236906;

	rk4_fin_weight = Array<double>{3};

	rk4_fin_weight(0) = 0.517231671970585; // u2
	rk4_fin_weight(1) = 0.096059710526147; // u3
	rk4_fin_weight(2) = 0.063692468666290; // du3

	ComputeTimeWeight();
}
*/
