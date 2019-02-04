/**
 * @file 	simulationarchive.h
 * @brief 	Tools for creating and readin a Simulation Archive binary file.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2016 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef SIMULATIONARCHIVE_H
#define SIMULATIONARCHIVE_H

#include <stdint.h>

struct reb_simulation;
struct reb_particles;

void reb_simulationarchive_heartbeat(struct reb_simulation* const r);  ///< Internal function to handle outputs for the Simulation Archive.
void reb_read_simulationarchive_with_messages(struct reb_simulationarchive* sa, const char* filename, enum reb_input_binary_messages* warnings); ///< Internal function to read one snapshot from a simulation archive.


#endif 	// SIMULATIONARCHIVE_H
