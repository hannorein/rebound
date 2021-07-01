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

#ifndef _UVARS_H_
#define _UVARS_H_
#include "Simulation.h"

// Build switch for using long double implementation of universal variables.
#ifdef USE_LONG_DOUBLES
  #define datatype_t long double
#else
  #define datatype_t double
#endif

typedef struct _StumpfCoefficients
{
  datatype_t * __restrict__ c0;
  datatype_t * __restrict__ c1;
  datatype_t * __restrict__ c2;
  datatype_t * __restrict__ c3;
}StumpfCoefficients;

typedef struct UNIVERSAL_VARS
{
  double * __restrict__ t0;
  double * __restrict__ tLast;
  double * uv_csq;
  double * uv_csp;
  double * uv_csv;  
  datatype_t * dt;
  datatype_t * __restrict__ Q0;
  datatype_t * __restrict__ V0;
  datatype_t * __restrict__ P0;
  datatype_t * __restrict__ Q1;
  datatype_t * __restrict__ V1;
  datatype_t * __restrict__ P1;  
  datatype_t * __restrict__ X;
  datatype_t * __restrict__ Q0_norm;
  datatype_t * __restrict__ beta;
  datatype_t * __restrict__ eta;
  datatype_t * __restrict__ zeta;
  datatype_t * __restrict__ period;
  datatype_t * __restrict__ Xperiod;
  uint32_t stateVectorSize;
  uint32_t stateVectorSize_l;
  uint32_t controlVectorSize;
  uint32_t controlVectorSize_l;

  StumpfCoefficients C;

  datatype_t mu;  /// G*mCentral
  datatype_t * mass; // Extra storage but means we dont have to cast everywhere.

  // Variables for storing classical orbital elements.
  double * e;
  double * a;
  double * h;
  double * h_norm;
  double * peri;
  double * apo;  
}UNIVERSAL_VARS;

//Interface functions.
void UniversalVars_Init(SIMULATION * z_sim);
void UniversalVars_Free(void);
void CalculateClassicalOrbitalElements(void);
void RebasisOsculatingOrbits_Momenta(double * z_Q, double * z_P, double z_t, uint32_t i);
void CalculateOsculatingOrbitsForSingleStep(double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis);
void ApplyCorrectorToOsculatingOrbitCalculation(double **Xosc_map, double t, uint32_t z_stagesPerStep);                                            


#endif