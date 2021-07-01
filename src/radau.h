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

#include "Simulation.h"
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
  double h;
  double t;

  // State vectors
  double * dX;
  double * dXtemp;
  double * dX0;
  double * dQ;
  double * dP;
  double * X;
  double * Q;
  double * P;

  // Stepsize control variables
  double aTol;
  double rTol;

  // Temporary buffers for rectifying before we output to file.
  double * Xout;
  double * Qout;
  double * Pout;

  // Step function pointers
  void (*step)(uint32_t *, double, double, uint32_t);
  void (*AnalyticalContinuation)(controlVars *, controlVars *, const double, const double, const uint32_t * const, const uint32_t);
  void (*CalculateGfromB)(void);
  double (*ReturnStepError)(double h, double t);
  void (*RejectStep)(void);

  uint32_t * rectifiedArray;

  // Buffer for predictors
  double * predictors;

  // Derivatives at start of the step.
  double * __restrict__ dState0;
  double * __restrict__ ddState0;

  double * dQ_dot_from_expansion;

  // Intermidiate derivatives.
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
  double * cs_state;
  double * cs_state_1st;
  double * csx;
  double * csv;

  double * csx_recti;
  double * csv_recti;

  double * cs_dq_dot;
  double * cs_dX;
  double * cs_dq;
  double * cs_dp;

  // Access pointers for testing
  controlVars * B;
  controlVars * Blast;

  controlVars * B_1st;
  controlVars * Blast_1st;
  controlVars * B_full;
  controlVars * Blast_full;  
  controlVars cs_B_full;

  // Variables for performance metrics
  uint64_t fCalls;
  uint64_t rectifications;
  uint64_t stepsTaken;
  uint32_t convergenceIterations;
  uint64_t stepsRejected;

  // File output fields
  double nextOutputTime;
  FILE * outputFile;
  double base;  // base value for the log scaling
  uint32_t output_samples_count;
  uint32_t t0_lim;


  double * b6_store;
  double * Xsize;

  // Terminate early function pointer.
  uint32_t (*fTerminate)(void);

  clock_t tStart;
  double cpuTimeUsed;

  double b6Max;
  double accMax;
  double * acc_ptr;

  double * q;
  double * p;
  double * q_dot;
  double * q_ddot;
  double * p_dot;

  double * q0;
  double * p0;
  double * q_dot0;
  double * q_ddot0;
  double * p_dot0;
}RADAU;

uint32_t Radau_integrate(void);
void Radau_Init(SIMULATION * z_sim);
void Radau_Free(void);
double Radau_CalculateStepSize(double h, double hLast, uint32_t * rejected, double t);


#endif
