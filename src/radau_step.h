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

#ifndef RADAU_STEP_H
#define RADAU_STEP_H
#include "Simulation.h"

void RadauStep15_Init(SIMULATION * z_sim);
void RadauStep15_Free(void);
void ControlVars_Free(controlVars * var);
void ControlVars_Init(controlVars * var, uint32_t size);
void ControlVars_Clear(controlVars * var);
void RejectStep(void);

void CalculatePredictors(double h, double hSample, double const * __restrict__ z_state0, double const * __restrict__ z_dState, 
                         double const * __restrict__ z_ddState, controlVars const * z_B, double * __restrict__ z_predictors, 
                         double const * __restrict__ z_csState, uint32_t const z_start, uint32_t const z_end, uint32_t const add_cs_var);
void CalculatePredictors_1stOrder(double h, double hSample, double * z_state0, double * z_dState, 
                                  controlVars * z_B, double * z_predictors, double *, uint32_t start, uint32_t length, uint32_t add_cs_var);

void CalculateNewState(double h, double * z_dState, double * z_ddState, controlVars * z_B, double * z_state0, double * z_csState, uint32_t z_start, uint32_t z_end);
void CalculateNewState_1stOrder(double h, double * z_dState, controlVars * z_B, double * z_state0, double * z_csState, uint32_t z_start, uint32_t z_end);
void AnalyticalContinuation(controlVars * z_B, controlVars * z_Blast, const double h, const double h_new, const uint32_t * const rectificationArray, const uint32_t step);
void CalculateGfromB(void);
#endif // RADAU_STEP_H
