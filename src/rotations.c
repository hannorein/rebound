/**
 * @file 	rotations.c
 * @brief 	Tools for rotations and quaternions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details This code uses the same conventions as the Apple SIMD quaternion framework. See e.g.: 
 *          https://github.com/xybp888/iOS-SDKs/blob/master/iPhoneOS13.0.sdk/usr/include/simd/quaternion.h
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

double reb_vec3d_dot(struct reb_vec3d a, struct reb_vec3d b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

double reb_vec3d_length_squared(struct reb_vec3d v){
    return reb_vec3d_dot(v, v);
}

struct reb_vec3d reb_vec3d_normalize(struct reb_vec3d v){
    return reb_vec3d_mul(v, 1./sqrt(reb_vec3d_length_squared(v)));
}

struct reb_quat reb_quat_mul(struct reb_quat p, struct reb_quat q){
    // v_rot = p * ( q * v)
    struct reb_quat r = {
        .r  = p.r*q.r  - p.ix*q.ix - p.iy*q.iy - p.iz*q.iz,
        .ix = p.r*q.ix + p.ix*q.r  + p.iy*q.iz - p.iz*q.iy,
        .iy = p.r*q.iy - p.ix*q.iz + p.iy*q.r  + p.iz*q.ix,
        .iz = p.r*q.iz + p.ix*q.iy - p.iy*q.ix + p.iz*q.r
    };
    return r;
}

double reb_quat_length_squared(struct reb_quat q){
    return q.r*q.r + q.ix*q.ix + q.iy*q.iy + q.iz*q.iz;
}
struct reb_quat reb_quat_identity(){
    struct reb_quat q = {.ix = 0.0, .iy = 0.0, .iz = 0.0, .r = 1.0 };
    return q;
}

struct reb_quat reb_quat_conjugate(struct reb_quat q){
    struct reb_quat c = {.ix = -q.ix, .iy = -q.iy, .iz = -q.iz, .r = q.r };
    return c;
}

struct reb_quat reb_quat_inverse(struct reb_quat q){
    struct reb_quat c = reb_quat_conjugate(q);
    double rl2 = 1./reb_quat_length_squared(q);
    c.r *= rl2;
    c.ix *= rl2;
    c.iy *= rl2;
    c.iz *= rl2;
    return c;
}

struct reb_quat reb_quat_init_with_angle_axis(double angle, struct reb_vec3d axis){
    axis = reb_vec3d_normalize(axis);
    double cos2 = cos(angle/2.0);
    double sin2 = sin(angle/2.0);
    struct reb_vec3d imag = reb_vec3d_mul(axis, sin2);
    struct reb_quat q = {.ix = imag.x, .iy = imag.y, .iz = imag.z, .r = cos2 };
    return q;
}

struct reb_quat reb_quat_init_with_orbital(const double Omega, const double inc, const double omega){
    // Murray and Dermot Eq. 2.121
    struct reb_vec3d x = {.x=1.0, .y=0.0, .z=0.0};
    struct reb_vec3d z = {.x=0.0, .y=0.0, .z=1.0};
    struct reb_quat P1 = reb_quat_init_with_angle_axis(omega, z);
    struct reb_quat P2 = reb_quat_init_with_angle_axis(inc, x);
    struct reb_quat P3 = reb_quat_init_with_angle_axis(Omega, z);
    return reb_quat_mul(P3, reb_quat_mul(P2, P1));
}

struct reb_vec3d reb_vec3d_rotate(struct reb_vec3d v, struct reb_quat q){
    struct reb_vec3d imag = reb_quat_imag(q);
    struct reb_vec3d t = reb_vec3d_mul(reb_vec3d_cross(imag,v), 2);
    return reb_vec3d_add(v, reb_vec3d_add(reb_vec3d_mul(t, q.r), reb_vec3d_cross(imag, t)));
}

static void reb_particle_irotate(struct reb_particle* p, struct reb_quat q){
    struct reb_vec3d pos = {p->x, p->y, p->z};
    pos = reb_vec3d_rotate(pos, q);
    p->x = pos.x;
    p->y = pos.y;
    p->z = pos.z;
    struct reb_vec3d vel = {p->vx, p->vy, p->vz};
    vel = reb_vec3d_rotate(vel, q);
    p->vx = vel.x;
    p->vy = vel.y;
    p->vz = vel.z;
}

void reb_simulation_irotate(struct reb_simulation* const sim, struct reb_quat q){
    const int N = sim->N;
    for (int i = 0; i < N; i++){
        struct reb_particle* p = &sim->particles[i];
        reb_particle_irotate(p,q);
    }
}


static inline struct reb_quat reb_quat_from_to_reduced(struct reb_vec3d from, struct reb_vec3d to) {
    // Internal use only
    struct reb_vec3d half = {.x=from.x+to.x, .y=from.y+to.y, .z=from.z+to.z};
    half = reb_vec3d_normalize(half);
    struct reb_vec3d cross = reb_vec3d_cross(from, half);
    double dot = reb_vec3d_dot(from, half);
    struct reb_quat q = {.ix=cross.x, .iy=cross.y, .iz=cross.z, .r=dot};
    return q;
}

struct reb_quat reb_quat_from_to(struct reb_vec3d from, struct reb_vec3d to) {

    if (reb_vec3d_dot(from, to) >= 0) {  // small angle
        return reb_quat_from_to_reduced(from, to);
    }

    //  More than 90 degrees apart, do rotation in two stages:
    //  (from -> half), (half -> to) 
    struct reb_vec3d half = {.x=from.x+to.x, .y=from.y+to.y, .z=from.z+to.z};
    half = reb_vec3d_normalize(half);

    if (reb_vec3d_length_squared(half) == 0) {
        //  half is nearly zero, so from and to point in nearly opposite directions
        //  and the rotation is numerically underspecified. Pick an axis orthogonal
        //  to the vectors, and use an angle of pi radians.
        struct reb_vec3d abs_from = {.x=fabs(from.x), .y=fabs(from.y), .z=fabs(from.z)};
        if (abs_from.x <= abs_from.y && abs_from.x <= abs_from.z){
            struct reb_vec3d axis = {.x=1, .y=0, .z = 0};
            axis = reb_vec3d_cross(from, axis);
            struct reb_quat q = {.ix=axis.x, .iy=axis.y, .iz=axis.z, .r=0.0};
            return q;
        }
        if (abs_from.y <= abs_from.z){
            struct reb_vec3d axis = {.x=0, .y=1, .z = 0};
            axis = reb_vec3d_cross(from, axis);
            struct reb_quat q = {.ix=axis.x, .iy=axis.y, .iz=axis.z, .r=0.0};
            return q;
        }
        struct reb_vec3d axis = {.x=0, .y=0, .z = 1};
        axis = reb_vec3d_cross(from, axis);
        struct reb_quat q = {.ix=axis.x, .iy=axis.y, .iz=axis.z, .r=0.0};
        return q;
    }

    return reb_quat_mul(reb_quat_from_to_reduced(from, half), reb_quat_from_to_reduced(half, to));
}

#define MIN_INC 1.e-8 
void reb_quat_to_orbital(struct reb_quat q, double* Omega, double* inc, double* omega){
    // see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0276302
    // and https://github.com/evbernardes/quaternion_to_euler/blob/main/euler_from_quat.py
    // Works but angles doen't always land in the right quadrant.
    double ap = q.r;
    double bp = q.iz;
    double cp = q.ix;
    double dp = q.iy;
    *inc = acos(2.0*(ap*ap+bp*bp) - 1.0);
    int safe1 =  (fabs(*inc) > MIN_INC);
    int safe2 =  (fabs(*inc - M_PI) > MIN_INC);

    if (safe1 && safe2){
        double half_sum = atan2(bp, ap);
        double half_diff = atan2(dp, cp);
        *omega = half_sum - half_diff;
        *Omega = half_sum + half_diff;
    }else{
        *Omega = 0;
        if (!safe1){
            double half_sum = atan2(bp, ap);
            *omega = 2.0 * half_sum;
        }else{
            double half_diff = atan2(dp, cp);
            *omega = 2.0 * half_diff;
        }
    }
    if (*omega < 0){
        *omega += M_PI*2.0;
    }
    if (*Omega < 0){
        *Omega += M_PI*2.0;
    }


}
