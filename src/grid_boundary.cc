/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "grid.hh"

void Grid::PeriodicBoundaryLeft()
{
	int i = 0;
	for (int m = 0; m < NQUANT; m++) {
		for (int k = 0; k < NGHOST; k++) {
			for (int j = NGHOST; j < nv-NGHOST; j++) {
				cons(m,i+k,j) = cons(m,nu-2*NGHOST+k,j);
			}
		}
	}
}

void Grid::PeriodicBoundaryRight()
{
	int i = nu-1;
	for (int m = 0; m < NQUANT; m++) {
		for (int k = 0; k < NGHOST; k++) {
			for (int j = NGHOST; j < nv-NGHOST; j++) {
				cons(m,i-k,j) = cons(m,2*NGHOST-1-k,j);
			}
		}
	}
}

void Grid::PeriodicBoundaryBot()
{
	int j = 0;
	for (int m = 0; m < NQUANT; m++) {
		for (int i = NGHOST; i < nu-NGHOST; i++) {
			for (int k = 0; k < NGHOST; k++) {
				cons(m,i,j+k) = cons(m,i,nv-2*NGHOST+k);
			}
		}
	}
}

void Grid::PeriodicBoundaryTop()
{
	int j = nv-1;
	for (int m = 0; m < NQUANT; m++) {
		for (int i = NGHOST; i < nu-NGHOST; i++) {
			for (int k = 0; k < NGHOST; k++) {
				cons(m,i,j-k) = cons(m,i,2*NGHOST-1-k);
			}
		}
	}
}

void Grid::PeriodicBoundaryLB()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,i,j) = cons(m,nu-2*NGHOST+i,nv-2*NGHOST+j);
			}
		}
	}
}

void Grid::PeriodicBoundaryRB()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,nu-NGHOST+i,j) = cons(m,i+NGHOST,nv-2*NGHOST+j);
			}
		}
	}
}

void Grid::PeriodicBoundaryRT()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,nu-NGHOST+i,nv-NGHOST+j) = cons(m,i+NGHOST,j+NGHOST);
			}
		}
	}
}

void Grid::PeriodicBoundaryLT()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,i,nv-NGHOST+j) = cons(m,nu-2*NGHOST+i,j+NGHOST);
			}
		}
	}
}

void Grid::InflowBoundaryLeft(double rho, double vx, double vy, double press)
{
	Array<double> tmp_prim{NQUANT};
	Array<double> tmp_cons{NQUANT};
	tmp_prim(0) = rho;
	tmp_prim(1) = vx;
	tmp_prim(2) = vy;
	tmp_prim(3) = press;
	for (int k = 0; k < NGHOST; k++) {
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			PointPrimToCons(tmp_prim, tmp_cons);
			cons(0,k,j) = tmp_cons(0);
			cons(1,k,j) = tmp_cons(1);
			cons(2,k,j) = tmp_cons(2);
			cons(3,k,j) = tmp_cons(3);
		}
	}
}

void Grid::InflowBoundaryRight(double rho, double vx, double vy, double press)
{
	int i = nu-1;
	Array<double> tmp_prim{NQUANT};
	Array<double> tmp_cons{NQUANT};
	tmp_prim(0) = rho;
	tmp_prim(1) = vx;
	tmp_prim(2) = vy;
	tmp_prim(3) = press;
	for (int k = 0;  k < NGHOST; k++) {
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			PointPrimToCons(tmp_prim, tmp_cons);
			prim(0,i-k,j) = tmp_cons(0);
			prim(1,i-k,j) = tmp_cons(1);
			prim(2,i-k,j) = tmp_cons(2);
			prim(3,i-k,j) = tmp_cons(3);
		}
	}
}

void Grid::InflowBoundaryBot(double rho, double vx, double vy, double press)
{
	Array<double> tmp_prim{NQUANT};
	Array<double> tmp_cons{NQUANT};
	tmp_prim(0) = rho;
	tmp_prim(1) = vx;
	tmp_prim(2) = vy;
	tmp_prim(3) = press;
	for (int i = NGHOST; i < nu-NGHOST; i++) {
		for (int k = 0;  k < NGHOST; k++) {
			PointPrimToCons(tmp_prim, tmp_cons);
			prim(0,i,k) = tmp_cons(0);
			prim(1,i,k) = tmp_cons(1);
			prim(2,i,k) = tmp_cons(2);
			prim(3,i,k) = tmp_cons(3);
		}
	}
}

void Grid::InflowBoundaryTop(double rho, double vx, double vy, double press)
{
	int j = nv-1;
	Array<double> tmp_prim{NQUANT};
	Array<double> tmp_cons{NQUANT};
	tmp_prim(0) = rho;
	tmp_prim(1) = vx;
	tmp_prim(2) = vy;
	tmp_prim(3) = press;
	for (int i = NGHOST; i < nu-NGHOST; i++) {
		for (int k = 0;  k < NGHOST; k++) {
			PointPrimToCons(tmp_prim, tmp_cons);
			prim(0,i,j-k) = tmp_cons(0);
			prim(1,i,j-k) = tmp_cons(1);
			prim(2,i,j-k) = tmp_cons(2);
			prim(3,i,j-k) = tmp_cons(3);
		}
	}
}

void Grid::SmoothBoundaryLeft()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int k = 0;  k < NGHOST; k++) {
			for (int j = NGHOST; j < nv-NGHOST; j++) {
				cons(m,k,j) = cons(m,NGHOST,j);
			}
		}
	}
}

void Grid::SmoothBoundaryRight()
{
	int i = nu-1;
	for (int m = 0; m < NQUANT; m++) {
		for (int k = 0;  k < NGHOST; k++) {
			for (int j = NGHOST; j < nv-NGHOST; j++) {
				cons(m,i-k,j) = cons(m,nu-NGHOST-1,j);
			}
		}
	}
}

void Grid::SmoothBoundaryBot()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = NGHOST; i < nu-NGHOST; i++) {
			for (int k = 0; k < NGHOST; k++) {
				cons(m,i,k) = cons(m,i,NGHOST);
			}
		}
	}
}

void Grid::SmoothBoundaryTop()
{
	int j = nv-1;
	for (int m = 0; m < NQUANT; m++) {
		for (int i = NGHOST; i < nu-NGHOST; i++) {
			for (int k = 0; k < NGHOST; k++) {
				cons(m,i,j-k) = cons(m,i,nv-NGHOST-1);
			}
		}
	}
}

void Grid::SmoothBoundaryLB()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,i,j) = cons(m,NGHOST,NGHOST);
			}
		}
	}
}

void Grid::SmoothBoundaryRB()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,nu-NGHOST+i,j) = cons(m,nu-1-NGHOST,NGHOST);
			}
		}
	}
}

void Grid::SmoothBoundaryRT()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,nu-NGHOST+i,nv-NGHOST+j) = cons(m,nu-1-NGHOST,nv-1-NGHOST);
			}
		}
	}
}

void Grid::SmoothBoundaryLT()
{
	for (int m = 0; m < NQUANT; m++) {
		for (int i = 0; i < NGHOST; i++) {
			for (int j = 0; j < NGHOST; j++) {
				cons(m,i,nv-NGHOST+j) = cons(m,NGHOST,nv-1-NGHOST);
			}
		}
	}
}

void Grid::ReflectingBoundaryLeft()
{
	for (int k = 0;  k < NGHOST; k++) {
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			cons(0,k,j) = cons(0,2*NGHOST-1-k,j);
			cons(1,k,j) = -1 * cons(1,2*NGHOST-1-k,j);
			cons(2,k,j) = cons(2,2*NGHOST-1-k,j);
			cons(3,k,j) = cons(3,2*NGHOST-1-k,j);
		}
	}
}

void Grid::ReflectingBoundaryRight()
{
	int i = nu-1;
	for (int k = 0;  k < NGHOST; k++) {
		for (int j = NGHOST; j < nv-NGHOST; j++) {
			cons(0,i-k,j) = cons(0,i-2*NGHOST+1+k,j);
			cons(1,i-k,j) = -1 * cons(1,i-2*NGHOST+1+k,j);
			cons(2,i-k,j) = cons(2,i-2*NGHOST+1+k,j);
			cons(3,i-k,j) = cons(3,i-2*NGHOST+1+k,j);
		}
	}
}

void Grid::ReflectingBoundaryBot()
{
	for (int i = NGHOST; i < nu-NGHOST; i++) {
		for (int k = 0;  k < NGHOST; k++) {
			cons(0,i,k) = cons(0,i,2*NGHOST-1-k);
			cons(1,i,k) = cons(1,i,2*NGHOST-1-k);
			cons(2,i,k) = -1 * cons(2,i,2*NGHOST-1-k);
			cons(3,i,k) = cons(3,i,2*NGHOST-1-k);
		}
	}
}

void Grid::ReflectingBoundaryTop()
{
	int j = nv-1;
	for (int i = NGHOST; i < nu-NGHOST; i++) {
		for (int k = 0;  k < NGHOST; k++) {
			cons(0,i,j-k) = cons(0,i,j-2*NGHOST+1+k);
			cons(1,i,j-k) = cons(1,i,j-2*NGHOST+1+k);
			cons(2,i,j-k) = -1 * cons(2,i,j-2*NGHOST+1+k);
			cons(3,i,j-k) = cons(3,i,j-2*NGHOST+1+k);
		}
	}
}

void Grid::ReflectingBoundaryLB()
{
	for (int i = 0; i < NGHOST; i++) {
		for (int j = 0; j < NGHOST; j++) {
			cons(0,i,j) = cons(0,2*NGHOST-1-i,2*NGHOST-1-j);
			cons(1,i,j) = -1 * cons(1,2*NGHOST-1-i,2*NGHOST-1-j);
			cons(2,i,j) = -1 * cons(2,2*NGHOST-1-i,2*NGHOST-1-j);
			cons(3,i,j) = cons(3,2*NGHOST-1-i,2*NGHOST-1-j);
		}
	}
}

void Grid::ReflectingBoundaryRB()
{
	for (int i = 0; i < NGHOST; i++) {
		for (int j = 0; j < NGHOST; j++) {
			cons(0,nu-1-i,j) = cons(0,nu-2*NGHOST+i,2*NGHOST-1-j);
			cons(1,nu-1-i,j) = -1 * cons(1,nu-2*NGHOST+i,2*NGHOST-1-j);
			cons(2,nu-1-i,j) = -1 * cons(2,nu-2*NGHOST+i,2*NGHOST-1-j);
			cons(3,nu-1-i,j) = cons(3,nu-2*NGHOST+i,2*NGHOST-1-j);
		}
	}
}

void Grid::ReflectingBoundaryRT()
{
	for (int i = 0; i < NGHOST; i++) {
		for (int j = 0; j < NGHOST; j++) {
			cons(0,nu-1-i,nv-1-j) = cons(0,nu-2*NGHOST+i,nv-2*NGHOST+j);
			cons(1,nu-1-i,nv-1-j) = -1 * cons(1,nu-2*NGHOST+i,nv-2*NGHOST+j);
			cons(2,nu-1-i,nv-1-j) = -1 * cons(2,nu-2*NGHOST+i,nv-2*NGHOST+j);
			cons(3,nu-1-i,nv-1-j) = cons(3,nu-2*NGHOST+i,nv-2*NGHOST+j);
		}
	}
}

void Grid::ReflectingBoundaryLT()
{
	for (int i = 0; i < NGHOST; i++) {
		for (int j = 0; j < NGHOST; j++) {
			cons(0,i,nv-1-j) = cons(0,2*NGHOST-1-i,nv-2*NGHOST+j);
			cons(1,i,nv-1-j) = -1 * cons(1,2*NGHOST-1-i,nv-2*NGHOST+j);
			cons(2,i,nv-1-j) = -1 * cons(2,2*NGHOST-1-i,nv-2*NGHOST+j);
			cons(3,i,nv-1-j) = cons(3,2*NGHOST-1-i,nv-2*NGHOST+j);
		}
	}
}
