/**
 * @file 	binarydiff.h
 * @brief 	Binary diff allows to compare binary snapshots.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2018 Hanno Rein
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
#ifndef _BINARYDIFF_H
#define _BINARYDIFF_H

// Compares two simulations, stores difference in buffer.
//
// output_option:
// - If set to 0, differences are written to bufp in the form of reb_binary_field structs. 
// - If set to 1, differences are printed on the screen. 
// - If set to 2, only the return value indicates any differences.
// - If set to 3, differences are written to bufp in a human readable form.
//
// returns value:  0 is returned if the simulations do not differ (are equal). 1 is return if they differ.

int reb_binary_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option);


#endif // _BINARYDIFF_H
