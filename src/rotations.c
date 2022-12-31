/**
 * @file 	rotations.c
 * @brief 	Tools for rotations and quaternions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2022 Hanno Rein
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
#include <stdlib.h>
#include <math.h>
#include "rebound.h"


struct reb_quat reb_quat_identity(){
    struct reb_quat q = {.ix = 0.0, .iy = 0.0, .iz = 0.0, .r = 1.0 };
    return q;
}

struct reb_vec3d reb_quat_imag(struct reb_quat q){
    struct reb_vec3d i = {
        .x = q.ix,
        .y = q.iy,
        .z = q.iz
    };
    return i;
}


struct reb_vec3d reb_vec3d_mul(struct reb_vec3d v, double s){
    struct reb_vec3d nv = {
        .x = s*v.x,
        .y = s*v.y,
        .z = s*v.z
    };
    return nv;
}

struct reb_vec3d reb_vec3d_add(struct reb_vec3d v, struct reb_vec3d w){
    struct reb_vec3d nv = {
        .x = v.x + w.x,
        .y = v.y + w.y,
        .z = v.z + w.z
    };
    return nv;
}

struct reb_vec3d reb_vec3d_cross(struct reb_vec3d a, struct reb_vec3d b){
    struct reb_vec3d c = {
        .x = a.y*b.z - a.z*b.y,
        .y = a.z*b.x - a.x*b.z,
        .z = a.x*b.y - a.y*b.x,
    };
    return c;
}

struct reb_vec3d reb_quat_act(struct reb_quat q, struct reb_vec3d v){
    struct reb_vec3d imag = reb_quat_imag(q);
    struct reb_vec3d t = reb_vec3d_mul(reb_vec3d_cross(imag,v), 2);
    return reb_vec3d_add(v, reb_vec3d_add(reb_vec3d_mul(t, q.r), reb_vec3d_cross(imag, t)));
}

