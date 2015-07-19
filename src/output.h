/**
 * @file 	output.h
 * @brief 	Output routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	If MPI is enabled, most functions output one file per
 * node. They automatically add a subscript _0, _1, .. to each file.
 * The user has to join them together manually. One exception is 
 * reb_output_append_velocity_dispersion() which only outputs one file.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#ifndef _OUTPUT_H
#define _OUTPUT_H
struct reb_simulation;

/**
 * Appends the velocity dispersion of the particles to an ASCII file.
 * @param filename Output filename.
 */
void reb_output_append_velocity_dispersion(struct reb_simulation* r, char* filename);

#ifdef PROFILING
/**
 * Profiling categories
 */
enum profiling_categories {
	PROFILING_CAT_INTEGRATOR,
	PROFILING_CAT_BOUNDARY,
	PROFILING_CAT_GRAVITY,
	PROFILING_CAT_COLLISION,
#ifdef OPENGL
	PROFILING_CAT_VISUALIZATION,
#endif // OPENGL
	PROFILING_CAT_NUM,
};
void profiling_start(void);
void profiling_stop(int cat);
#define PROFILING_START() profiling_start();
#define PROFILING_STOP(C) profiling_stop(C);
#else // PROFILING
#define PROFILING_START()	// Dummy, do nothing 
#define PROFILING_STOP(C)	
#endif // PROFILING

#endif
