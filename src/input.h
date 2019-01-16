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

void reb_read_dp7(struct reb_dp7* dp7, const int N3, FILE* inf, char **restrict mem_stream); ///< Internal function to read dp7 structs from file.
int reb_input_field(struct reb_simulation* r, FILE* inf, enum reb_input_binary_messages* warnings, char **restrict mem_stream); ///< Read one field from inf stream into r. 

struct reb_simulation* reb_input_process_warnings(struct reb_simulation* r, enum reb_input_binary_messages warnings); ///< Process warning messages and print them on screen.

#define _INPUT_H

#endif
