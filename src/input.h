/**
 * @file 	input.h
 * @brief 	Parse command line options and read retart files.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#ifndef _INPUT_H
#define _INPUT_H

/**
 * Checks if a restart file was given as a command line argument. 
 * @return Returns 1 if simulation was restarted. 0 otherwise.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 */ 
int input_check_restart(int argc, char** argv);

/**
 * Reads a binary file.
 * @param filename Filename to be read.
 */
void input_binary(char* filename);

/**
 * Reads arguments from the command line.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @param argument Argument to look for.
 * @return Returns NULL if argument was not given. Return the argument otherwise.
 */
char* input_get_argument(int argc, char** argv, const char* argument);


#endif
