/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <thread>
#include <memory>
#include <vector>
#include <pthread.h>

#include "grid.h"
#include "riemann.h"
#include "integrator.h"
#include "broadcast.h"

double global_time;
double dt;
double step_time;
double step_dt;

double out_time;

class IntegratorThread {
public:
	int tid;
	Integrator &integrator;
	std::unique_ptr<std::thread> thread;
	Grid &global_grid;
	Grid local_grid;
	pthread_barrier_t *barrier;

	IntegratorThread(int tid, Integrator &integrator, Grid &g, pthread_barrier_t *barrier)
	: tid{tid}, integrator{integrator}, global_grid{g}, local_grid{global_time, dt, step_time, step_dt} {
		this->barrier = barrier;

		int ni_per_thread = (global_grid.nu - 2*NGHOST + NTHREAD - 1) / NTHREAD;

		local_grid.il = NGHOST + ni_per_thread * tid;
		if (tid == 0) {
			local_grid.ilr = local_grid.il - 1;
		} else {
			local_grid.ilr = local_grid.il;
		}

		if (tid == NTHREAD - 1) {
			local_grid.iu = global_grid.nu - NGHOST;
			local_grid.iuf = local_grid.iu + 1;
			local_grid.iur = local_grid.iu + 1;
		} else {
			local_grid.iu = NGHOST + ni_per_thread * (tid + 1);
			local_grid.iuf = local_grid.iu;
			local_grid.iur = local_grid.iu;
		}

		local_grid.jl = NGHOST;
		local_grid.ju = global_grid.nv - NGHOST;

		local_grid.tid = tid;
		local_grid.AttachReference(global_grid);
		start();
	}
	~IntegratorThread() {
		join();
		local_grid.DetachReference();
	}

	void start() {
		thread = std::make_unique<std::thread>(&IntegratorThread::thread_main, this);
	}
	void join() {
		if (thread) {
			if (thread->joinable()) {
				thread->join();
			}
			thread.reset();
		}
	}

	void take_timestep() {
		Array<double> *J;
		int &s = integrator.s;
		if (tid == 0) {
			s = 0;
			global_grid.cons_gen.copy_data_from(global_grid.cons);
			dt = DBL_MAX;
		}
		pthread_barrier_wait(barrier);

		while (s < integrator.nstep) {
			for (int dir = 0; dir < 2; dir++) {
				if (dir == 0) {
					J = &local_grid.Ju;
				} else {
					J = &local_grid.Jv;
				}

				local_grid.Reconstruct(dir);
				pthread_barrier_wait(barrier);

				local_grid.PrimLim(local_grid.Lprim);
				local_grid.PrimLim(local_grid.Rprim);
				local_grid.PrimToCons(local_grid.Lprim, local_grid.Lcons);
				local_grid.PrimToCons(local_grid.Rprim, local_grid.Rcons);
				pthread_barrier_wait(barrier);

				local_grid.Wavespeed(dir);

				// timestep determination
				pthread_barrier_wait(barrier);
				if (tid == 0 && s == 0) {
					global_grid.DetermineDt(dir);
				}
				pthread_barrier_wait(barrier);

				riemann::HLLC(local_grid.Lprim, local_grid.Lcons,
				local_grid.Lw, local_grid.Rprim, local_grid.Rcons, local_grid.Rw, *J, dir,
				local_grid.il, local_grid.iuf, local_grid.jl, local_grid.ju);
			}

			// finalize timestep determination
			if (tid == 0) {
				if (s == 0) {
					dt *= integrator.cfl_num;
				}
				step_dt = integrator.time_weight(s) * dt;
				if (s == 0) {
					step_time = global_time;
				} else {
					step_time = global_time + integrator.time_weight(s-1) * dt;
				}
			}
			pthread_barrier_wait(barrier);

			// hydro
			local_grid.CalculateFluxDiv();
			local_grid.CalculateSrc();

			integrator.AddFluxDivSrc(&local_grid);
			pthread_barrier_wait(barrier);

			local_grid.ConsLim();
			local_grid.ConsToPrim();

			pthread_barrier_wait(barrier);
			if (tid == 0) {
				global_grid.Boundary(step_time);
			}
			pthread_barrier_wait(barrier);

			local_grid.ConsLim();
			local_grid.ConsToPrim();

			if (tid == 0) {
				s++;
			}
			pthread_barrier_wait(barrier);
		}

		if (tid == 0) {
			global_time += dt;
		}
		pthread_barrier_wait(barrier);
	}

	void thread_main() {
		for (int epoch = 0; epoch < integrator.max_epoch; epoch++) {
			if (tid == 0) {
				printf("t = %.3e\tdt = %.3e\t%.2f%%\n", global_time, dt, 100*global_time/integrator.out_tf);
				if (global_time >= out_time) {
					global_grid.mutex.lock();
					global_grid.broadcast_signal = true;
					global_grid.mutex.unlock();
					global_grid.cond.notify_all();
					out_time = global_time + integrator.out_dt;
				}
			}

			take_timestep();
		}
	}
};

int main()
{
	std::vector<std::unique_ptr<IntegratorThread>> integrator_threads;
	pthread_barrier_t barrier;

	Grid global_grid{global_time, dt, step_time, step_dt};
	global_grid.tid = -1;
	global_grid.InitGrid();

	Integrator integrator;
	integrator.Property();

	Broadcast broadcast{global_grid, 9743, 2, 0, 24};

	pthread_barrier_init(&barrier, NULL, NTHREAD);
	for (int tid = 0; tid < NTHREAD; tid++) {
		integrator_threads.push_back(std::make_unique<IntegratorThread>(tid, integrator, global_grid, &barrier));
	}

	for (int tid = 0; tid < NTHREAD; tid++) {
		integrator_threads[tid]->join();
	}
	broadcast.join();

	pthread_barrier_destroy(&barrier);
	return 0;
}
