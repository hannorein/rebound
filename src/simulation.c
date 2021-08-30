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

#include <stdlib.h>
#include <string.h>
#include "rebound.h"

static SIMULATION * sim = NULL;

SIMULATION * Simulation_Init(struct reb_simulation* r, uint32_t z_n)
{
  sim = (SIMULATION *)malloc(sizeof(SIMULATION));

  // Initialise all object variables to zero to begin with.
  memset(sim, 0, sizeof(SIMULATION));

  // Set control variables to initial values
  r->ri_tes.stateVectorLength = 2*3*r->N;
  r->ri_tes.stateVectorSize = r->ri_tes.stateVectorLength * sizeof(double);
  r->ri_tes.controlVectorSize = r->N * sizeof(double);

  // Allocate memory
  sim->mass = (double *)malloc(r->ri_tes.controlVectorSize);
  sim->X_dh = (double *)malloc(r->ri_tes.stateVectorSize);
  sim->Q_dh = sim->X_dh;
  sim->P_dh = &sim->X_dh[r->ri_tes.stateVectorLength/2];

  // Ensure we are clean for each integration.
  memset(sim->mass, 0, r->ri_tes.controlVectorSize);
  memset(sim->X_dh, 0, r->ri_tes.stateVectorSize);

  sim->termination_check_enable = 0;

  return sim;
}


void Simulation_Free(void)
{
  free(sim->X_dh);
  free(sim->mass);
  free(sim);
}
