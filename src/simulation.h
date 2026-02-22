/**
 * @file 	simulation.h
 * @brief 	Functions which create, free, and operate on simulations.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2026 Hanno Rein
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
#ifndef SIMULATION_H
#define SIMULATION_H

// Finds the two largest particles in the simulation. *p1 and *p2 will be set to the indices of the largest particles.
void reb_simulation_two_largest_particles(struct reb_simulation* r, int* p1, int* p2);

// Only free memory in pointers of a simulation, but not the simulation itself.
DLLEXPORT void reb_simulation_free_pointers(struct reb_simulation* const r);

// Used internally and by python. Should not be called by the user.
DLLEXPORT void reb_simulation_init(struct reb_simulation* r); 

// Used from python and for display.
DLLEXPORT void reb_simulation_copy_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum reb_simulation_binary_error_codes* warnings);

// Serialization functions. Caller is responsible for allocating memory. Null pointers will not be set/read. 
DLLEXPORT void reb_simulation_get_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]);
DLLEXPORT void reb_simulation_set_serialized_particle_data(struct reb_simulation* r, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]);

// Used for unit testing only.
DLLEXPORT size_t reb_simulation_struct_size();

// Python workaround for function pointer.
DLLEXPORT void reb_simulation_set_collision_resolve(struct reb_simulation* r, enum REB_COLLISION_RESOLVE_OUTCOME (*resolve) (struct reb_simulation* const r, struct reb_collision c));

// Return the difference between two simulations as a human readable string. Returned pointer needs to be freed by caller.
DLLEXPORT char* reb_simulation_diff_char(struct reb_simulation* r1, struct reb_simulation* r2);

#endif 	// SIMULATION_H
