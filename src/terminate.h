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

#ifndef TERMINATE_H
#define TERMINATE_H

#include "Simulation.h"

typedef enum
{
  none = 0,
  intersection = 1,
  collision = 2,
  combined = 3
}TERMINATION;

void Terminate_Init(SIMULATION * sim, uint32_t (**f)(void), TERMINATION);
void Terminate_Free(SIMULATION * z_sim);

#endif
