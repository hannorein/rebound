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

#include "rebound.h"

typedef struct _StumpfCoefficients
{
  double * __restrict__ c0;
  double * __restrict__ c1;
  double * __restrict__ c2;
  double * __restrict__ c3;
}StumpfCoefficients;

typedef struct UNIVERSAL_VARS
{
  double * __restrict__ t0;
  double * __restrict__ tLast;
  double * uv_csq;
  double * uv_csp;
  double * uv_csv;  
  double * dt;
  double * __restrict__ Q0;
  double * __restrict__ V0;
  double * __restrict__ P0;
  double * __restrict__ Q1;
  double * __restrict__ V1;
  double * __restrict__ P1;  
  double * __restrict__ X;
  double * __restrict__ Q0_norm;
  double * __restrict__ beta;
  double * __restrict__ eta;
  double * __restrict__ zeta;
  double * __restrict__ period;
  double * __restrict__ Xperiod;
  uint32_t stateVectorSize;
  uint32_t controlVectorSize;

  StumpfCoefficients C;

  double mu;  /// G*mCentral

  // Variables for storing classical orbital elements.
  double * e;
  double * a;
  double * h;
  double * h_norm;
  double * peri;
  double * apo;  
}UNIVERSAL_VARS;

//Interface functions.
void UniversalVars_Init(struct reb_simulation* const r);
void UniversalVars_Free(struct reb_simulation* r);
void CalculateClassicalOrbitalElements(struct reb_simulation* r);
void RebasisOsculatingOrbits_Momenta(struct reb_simulation* r, double * z_Q, double * z_P, double z_t, uint32_t i);
void CalculateOsculatingOrbitsForSingleStep(struct reb_simulation* r, double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis);
void ApplyCorrectorToOsculatingOrbitCalculation(struct reb_simulation* r, double **Xosc_map, double t, uint32_t z_stagesPerStep);                                            


#endif