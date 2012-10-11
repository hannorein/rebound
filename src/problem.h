/**
 * @file 	problem.h
 * @brief 	The user must implement all functions declared in this file. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This file contains the function declarations that must be 
 * implemented by the user in the problem.c file. All functions must be 
 * present, even if they are empty and have no use.
 * 
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
#ifndef _PROBLEM_H
#define _PROBLEM_H
/**
 * Main initialization function.
 * In this function, the user should set the timestep, the boxsize and
 * set up all particles according to the initial conditions and add them
 * to the simulation using the particles_add() routine. The user may also 
 * read in command line arguments to change parameters.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 */ 
void problem_init(int argc, char* argv[]);
/**
 * This function is called at the beginning of the simulation, at the end of
 * each timestep and at the end of the simulation. The user can call either
 * generic output routines such as output_ascii() or create their own problem
 * specific output routines.
 */
void problem_output();
/**
 * This function is called in the middle of each timestep (after the gravitational
 * acceleration is added to the particles. The user may add 
 * any problem specific work (additional forces, etc) here.
 */
void problem_inloop();
/**
 * This function is called at the end of the simulation when t>=tmax.
 * Note that it is not called when the simulation stopped for another 
 * reason (e.g. user interaction or crash). 
 */ 
void problem_finish();

/*
 * This function allows the user to add additional (non-gravitational) forces.
 */
extern void (*problem_additional_forces) ();
#endif //_PROBLEM_H
