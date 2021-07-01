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

#ifndef _DHEM_H_
#define _DHEM_H_

#include "Simulation.h"

typedef enum _RHS_CONFIG_
{
  rhs_full = 0,         // The standard DHEM RHS.
  rhs_extended = 1      // The DHEM RHS but with extended precision for the osculating orbit.
} RHS_CONFIG;

typedef struct DHEM
{
  double * Xosc;
  double * Qosc;
  double * Posc;
  double * Vosc;

  double * XoscStore;
  double ** XoscArr;

  double * XoscPredStore;
  double ** XoscPredArr;  

  // CS variables for the osculating orbits.
  double * Xosc_cs;
  double * XoscStore_cs;
  double ** XoscArr_cs;
  double * Qosc_cs;
  double * Posc_cs;

  double * Xosc_dotStore;
  double ** Xosc_dotArr;

  double * Xosc_dot;
  double * Qosc_dot;
  double * Posc_dot;

  double * X;
  double * Q;
  double * P;

  double * X_dot;
  double * Q_dot;
  double * P_dot;

  double * rectifyTimeArray;  /// The time at which we need to rectify each body.

  double * __restrict__ m;
  double * __restrict__ m_inv;
  double mTotal;
  double n;

  double * rectificationPeriod;

  double * dQdot_cs;
  double * dQddot_cs;
  double * dPdot_cs;

  // Long double implementation storage for osculating orbits.
  long double * XoscStore_ld;
  long double ** XoscArr_ld;
  long double * Xosc_dotStore_ld;
  long double ** Xosc_dotArr_ld;  


  long double * Xosc_ld;
  long double * Qosc_ld;
  long double * Posc_ld;
  long double * Vosc_ld;
  long double * Xosc_dot_ld;
  long double * Qosc_dot_ld;
  long double * Posc_dot_ld;

  uint32_t final_stage_index;

  RHS_CONFIG rhsConfig;  
}DHEM;

void dhem_CalcOscOrbitsForAllStages(double t0, double h, double * hArr, uint32_t z_stagesPerStep, uint32_t z_rebasis);
double dhem_CalculateHamiltonian(double * Q, double * P);
void dhem_ConvertToDHCoords(double * Q, double * V, double * Qout, double * Pout);
void dhem_ConvertToCOM(double * Q, double * V, double * Qout, double * Vout);
void dhem_rhs(double const * __restrict__ const dQ, double const * __restrict__ const dP, double * __restrict__ const dQ_dot,
              double * __restrict__ const dP_dot, double * __restrict__ const dQ_ddot, double * __restrict__ const dP_ddot);
void dhem_rhs_full(double * r, double * acc);
void dhem_rhs_wrapped(double * dQ, double * dP, double * dQ_dot,
                      double * dP_dot, double * dQ_ddot, double * dP_ddot, uint32_t stageNumber,
                      double * cs1, double * cs2);      
void dhem_CalculateOsculatingOrbitDerivatives_Momenta(double const * const __restrict__ Qosc, double const * const __restrict__ Posc, 
                                                      double * const __restrict__ Qosc_dot, double * const __restrict__ Posc_dot);
uint32_t dhem_RectifyOrbits(double t, double * Q, double * P,
                            double * dQ, double * dP, uint32_t * rectifiedArray, uint32_t stageNumber);
void dhem_Init(SIMULATION * z_sim, double z_rectificationPeriodDefault, uint32_t z_stagesPerStep);
void dhem_Free(void);
void dh_full_rhs(double * Q, double * P, double * Qdot, double * Qddot, double * Pdot, double * mass, uint32_t n);

#endif
