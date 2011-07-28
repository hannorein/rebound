/**
 * @file main.h
 * @author Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of nbody.
 *
 * nbody is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * nbody is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with nbody.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef _MAIN_H
#define _MAIN_H

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

extern double softening;
extern double G;
extern double t;
extern double tmax;
extern double dt;
extern int N;
extern int N_active_first;
extern int N_active_last;
extern double timing_initial;
extern int root_nx;
extern int root_ny;
extern int root_nz;
extern int root_n;
extern double boxsize;
extern double boxsize_x;
extern double boxsize_y;
extern double boxsize_z;
extern double boxsize_max;

void init_box();
void iterate();

#endif
