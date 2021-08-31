// TES is an open source integration package for modelling exoplanet evolution.
// Copyright (C) <2021>  <Peter Bartram, Alexander Wittig>

// TES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>. 

#ifndef _IAS15_H_
#define _IAS15_H_

#include "rebound.h"
#include <stdio.h>
#include <time.h>

typedef struct _controlVars {
    double* __restrict__ p0;
    double* __restrict__ p1;
    double* __restrict__ p2;
    double* __restrict__ p3;
    double* __restrict__ p4;
    double* __restrict__ p5;
    double* __restrict__ p6;
    uint32_t size;
}controlVars;


typedef struct RADAU
{
  // State vectors
  double * dX;
  double * dQ;
  double * dP;
  // Buffers for rectifying into before performing a synchronisation.
  double * Xout;
  double * Qout;
  double * Pout;
  uint32_t * rectifiedArray;
  // Buffer for predictors
  double * predictors;
  // Derivatives at start of the step.
  double * __restrict__ dState0;
  double * __restrict__ ddState0;
  // Intermediate derivatives.
  double * __restrict__ dState;
  double * __restrict__ ddState;
  // Compensated summation arrays for gravity
  double * __restrict__ cs_dState0;
  double * __restrict__ cs_ddState0;
  double * __restrict__ cs_dState;
  double * __restrict__ cs_ddState;
  // Compensated summation arrays for B's.
  controlVars cs_B;
  controlVars cs_B1st;
  // Compensated summation array for the predictor and corrector.
  double * cs_dX;
  double * cs_dq;
  double * cs_dp;
  // Integrator coefficients.
  controlVars * B;
  controlVars * Blast;
  controlVars * B_1st;
  controlVars * Blast_1st; 
  // Variables for performance metrics
  uint64_t fCalls;
  uint64_t rectifications;
  uint32_t convergenceIterations;
  // Iteration convergence variables.
  double * b6_store;
  double * acc_ptr;
}RADAU;

void Radau_Init(struct reb_simulation* r);
void Radau_Free(struct reb_simulation* r);
double Radau_CalculateStepSize(struct reb_simulation* r, double h, double hLast, double t);
void ClearRectifiedBFields(struct reb_simulation* r, controlVars * B, uint32_t * rectifiedArray);
double Radau_SingleStep(struct reb_simulation* r, double z_t, double dt, double dt_last_done);
#endif
