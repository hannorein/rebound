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

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <stdint.h>
#include <stdio.h>

typedef struct UNIVERSAL_VARS UNIVERSAL_VARS;
typedef struct DHEM DHEM;
typedef struct RADAU RADAU;

/**
 * This is the main error enum for the simulation.
 */
typedef enum ERROR_SIMULATION_
{
	SIM_ERROR_OK = 0,	/// Everything is OK.

}ERROR_SIMULATION;


typedef struct _SIMULATION_
{
	// ##############################################
	// To use rebound data structure instead
	// ##############################################
	double t0;								/// Initial time
	uint32_t n;								/// Number of particles presently
	double rTol;								/// Relative tolerance used by the integrator.
	double dQcutoff;						/// Delta Q growth allowed before rectification.

	// ##############################################
	// To keep / refactor
	// ##############################################

	uint32_t stateVectorLength;		/// Length of the state vector in doubles.
	uint32_t n3;					/// 3*number of bodies.
	uint32_t stateVectorSize;		/// Size in bytes of the state vector
	uint32_t stateVectorSize_ld;	/// Size of state vector in long double data type.
	uint32_t controlVectorSize; /// Size in bytes of n * sizeof(double)
	uint32_t controlVectorLength;		/// Length of the control vector in doubles.
	UNIVERSAL_VARS * uVars;		/// Pointer to the universal variables module
	DHEM * rhs;								/// Pointer to the DHEM rhs
	RADAU * radau;  /// Pointer to our integrator

	double * mass;							/// Initial particle masses
	double * X_dh;						/// Memory for current state in dh coords.
	double * Q_dh;						/// Current state in dh coords.
	double * P_dh;						/// Current state in dh coords.

	// ##############################################
	// To throw away
	// ##############################################
	uint32_t step;
	uint32_t rectificationCount;
	uint32_t h_last_done;
	uint32_t termination_check_enable;

	// RHS function pointers to be provided by all force models.
	void (*f_rhs)(struct reb_simulation* r, double * dQ, double * dP, double * dQ_dot,
					double * dP_dot, double * dQ_ddot, double * dP_ddot,
									uint32_t stageNumber, double * cs1, double * cs2);
	void (*f_rhs_full)(double * r, double * acc);								
	void (*fStartOfStep)(struct reb_simulation* r, double t0, double h, double * hArr, uint32_t z_stagesPerStep, uint32_t z_rebasis);
	uint32_t (*fRectify)(double t, double * Q, double * P,
								double * dQ, double * dP, uint32_t * rectifiedArray, uint32_t stageNumber);
		void (*fPerformSummation)(double *, double *, double *, double *, uint32_t);
		double (*fCalculateInvariant)(struct reb_simulation* r, double * Q, double * P);
	void (*fCalculateOsculatingOrbitNorm)(double * Xosc_norm);
	void (*fCalculate_dQdot)(double * dP, double * dQdot, double * Posc);
	void (*f_dh_full_rhs)(double * Q, double * P, double * Qdot, double * Qddot, double * Pdot, double * mass, uint32_t n);

	
}SIMULATION;


SIMULATION * Simulation_Init(uint32_t z_n);
void Simulation_Free(void);

#endif /* SIMULATION_H_ */
