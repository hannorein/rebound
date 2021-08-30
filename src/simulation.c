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

#include "simulation.h"
#include <stdlib.h>
#include <string.h>
#include "rebound.h"

static SIMULATION * sim = NULL;

SIMULATION * Simulation_Init(uint32_t z_n)
{
  sim = (SIMULATION *)malloc(sizeof(SIMULATION));

  // Initialise all object variables to zero to begin with.
  memset(sim, 0, sizeof(SIMULATION));

  // Set control variables to initial values
  sim->n = z_n;
  sim->n3 = 3*z_n;
  sim->stateVectorLength = 2*3*sim->n;
  sim->stateVectorSize = sim->stateVectorLength * sizeof(double);
  sim->stateVectorSize_ld = sim->stateVectorLength * sizeof(long double);
  sim->controlVectorSize = sim->n * sizeof(double);

  // Allocate memory
  sim->mass = (double *)malloc(sim->controlVectorSize);
  sim->X_dh = (double *)malloc(sim->stateVectorSize);
  sim->Q_dh = sim->X_dh;
  sim->P_dh = &sim->X_dh[sim->stateVectorLength/2];

  // Ensure we are clean for each integration.
  memset(sim->mass, 0, sim->controlVectorSize);
  memset(sim->X_dh, 0, sim->stateVectorSize);

  sim->termination_check_enable = 0;

  return sim;
}


void Simulation_Free(void)
{
  free(sim->X_dh);
  free(sim->mass);
  free(sim);
}
