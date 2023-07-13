/**
 * @file    transformations.c
 * @brief   Transformations back and forth between different coordinate systems.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file collects all the transformations used by the different integrators
 * between different coordinate systems.
 *
 * @section     LICENSE
 * Copyright (c) 2017 Hanno Rein, Dan Tamayo.
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

#include "transformations.h"
#include "rebound.h"

/******************************
 * Jacobi */

void reb_transformations_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta = p_mass[0].m;
    double s_x = eta * particles[0].x;
    double s_y = eta * particles[0].y;
    double s_z = eta * particles[0].z;
    double s_vx = eta * particles[0].vx;
    double s_vy = eta * particles[0].vy;
    double s_vz = eta * particles[0].vz;
    for (unsigned int i=1;i<N_active;i++){
        const double ei = 1./eta;
        const struct reb_particle pi = particles[i];
        eta += p_mass[i].m;
        const double pme = eta*ei;
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
        s_x  = s_x  * pme + p_mass[i].m*p_j[i].x ;
        s_y  = s_y  * pme + p_mass[i].m*p_j[i].y ;
        s_z  = s_z  * pme + p_mass[i].m*p_j[i].z ;
        s_vx = s_vx * pme + p_mass[i].m*p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m*p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m*p_j[i].vz;
    }
    const double ei = 1./eta;
    for (unsigned int i=N_active;i<N;i++){
        const struct reb_particle pi = particles[i];
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
    }
    const double Mtotal  = eta;
    const double Mtotali = 1./Mtotal;
    p_j[0].m = Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
}

void reb_transformations_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta = p_mass[0].m;
    double s_x = eta * particles[0].x;
    double s_y = eta * particles[0].y;
    double s_z = eta * particles[0].z;
    double s_vx = eta * particles[0].vx;
    double s_vy = eta * particles[0].vy;
    double s_vz = eta * particles[0].vz;
    double s_ax = eta * particles[0].ax;
    double s_ay = eta * particles[0].ay;
    double s_az = eta * particles[0].az;
    for (unsigned int i=1;i<N_active;i++){
        const double ei = 1./eta;
        const struct reb_particle pi = particles[i];
        eta += p_mass[i].m;
        const double pme = eta*ei;
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
        s_x  = s_x  * pme + p_mass[i].m*p_j[i].x ;
        s_y  = s_y  * pme + p_mass[i].m*p_j[i].y ;
        s_z  = s_z  * pme + p_mass[i].m*p_j[i].z ;
        s_vx = s_vx * pme + p_mass[i].m*p_j[i].vx;
        s_vy = s_vy * pme + p_mass[i].m*p_j[i].vy;
        s_vz = s_vz * pme + p_mass[i].m*p_j[i].vz;
        s_ax = s_ax * pme + p_mass[i].m*p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m*p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m*p_j[i].az;
    }
    const double ei = 1./eta;
    for (unsigned int i=N_active;i<N;i++){
        const struct reb_particle pi = particles[i];
        p_j[i].m = pi.m;
        p_j[i].x = pi.x - s_x*ei;
        p_j[i].y = pi.y - s_y*ei;
        p_j[i].z = pi.z - s_z*ei;
        p_j[i].vx = pi.vx - s_vx*ei;
        p_j[i].vy = pi.vy - s_vy*ei;
        p_j[i].vz = pi.vz - s_vz*ei;
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
    }
    const double Mtotal  = eta;
    const double Mtotali = 1./Mtotal;
    p_j[0].m = Mtotal;
    p_j[0].x = s_x * Mtotali;
    p_j[0].y = s_y * Mtotali;
    p_j[0].z = s_z * Mtotali;
    p_j[0].vx = s_vx * Mtotali;
    p_j[0].vy = s_vy * Mtotali;
    p_j[0].vz = s_vz * Mtotali;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

void reb_transformations_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta = p_mass[0].m;
    double s_ax = eta * particles[0].ax;
    double s_ay = eta * particles[0].ay;
    double s_az = eta * particles[0].az;
    for (unsigned int i=1;i<N_active;i++){
        const double ei = 1./eta;
        const struct reb_particle pi = particles[i];
        eta += p_mass[i].m;
        const double pme = eta*ei;
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
        s_ax = s_ax * pme + p_mass[i].m*p_j[i].ax;
        s_ay = s_ay * pme + p_mass[i].m*p_j[i].ay;
        s_az = s_az * pme + p_mass[i].m*p_j[i].az;
    }
    const double ei = 1./eta;
    for (unsigned int i=N_active;i<N;i++){
        const struct reb_particle pi = particles[i];
        p_j[i].ax = pi.ax - s_ax*ei;
        p_j[i].ay = pi.ay - s_ay*ei;
        p_j[i].az = pi.az - s_az*ei;
    }
    const double Mtotal  = eta;
    const double Mtotali = 1./Mtotal;
    p_j[0].ax = s_ax * Mtotali;
    p_j[0].ay = s_ay * Mtotali;
    p_j[0].az = s_az * Mtotali;
}

void reb_transformations_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta  = p_j[0].m;
    double s_x  = p_j[0].x  * eta;
    double s_y  = p_j[0].y  * eta;
    double s_z  = p_j[0].z  * eta;
    double s_vx = p_j[0].vx * eta;
    double s_vy = p_j[0].vy * eta;
    double s_vz = p_j[0].vz * eta;
    for (unsigned int i=N-1;i>=N_active;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        particles[i].x  = pji.x  + s_x  * ei;
        particles[i].y  = pji.y  + s_y  * ei;
        particles[i].z  = pji.z  + s_z  * ei;
        particles[i].vx = pji.vx + s_vx * ei;
        particles[i].vy = pji.vy + s_vy * ei;
        particles[i].vz = pji.vz + s_vz * ei;
    }
    for (unsigned int i=N_active-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        s_x  = (s_x  - p_mass[i].m * pji.x ) * ei;
        s_y  = (s_y  - p_mass[i].m * pji.y ) * ei;
        s_z  = (s_z  - p_mass[i].m * pji.z ) * ei;
        s_vx = (s_vx - p_mass[i].m * pji.vx) * ei;
        s_vy = (s_vy - p_mass[i].m * pji.vy) * ei;
        s_vz = (s_vz - p_mass[i].m * pji.vz) * ei;
        particles[i].x  = pji.x  + s_x ;
        particles[i].y  = pji.y  + s_y ;
        particles[i].z  = pji.z  + s_z ;
        particles[i].vx = pji.vx + s_vx;
        particles[i].vy = pji.vy + s_vy;
        particles[i].vz = pji.vz + s_vz;
        eta -= p_mass[i].m;
        s_x  *= eta;
        s_y  *= eta;
        s_z  *= eta;
        s_vx *= eta;
        s_vy *= eta;
        s_vz *= eta;
    }
    const double mi = 1./eta;
    particles[0].x  = s_x  * mi;
    particles[0].y  = s_y  * mi;
    particles[0].z  = s_z  * mi;
    particles[0].vx = s_vx * mi;
    particles[0].vy = s_vy * mi;
    particles[0].vz = s_vz * mi;
}

void reb_transformations_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta  = p_j[0].m;
    double s_x  = p_j[0].x  * eta;
    double s_y  = p_j[0].y  * eta;
    double s_z  = p_j[0].z  * eta;
    for (unsigned int i=N-1;i>=N_active;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        particles[i].x  = pji.x  + s_x*ei ;
        particles[i].y  = pji.y  + s_y*ei ;
        particles[i].z  = pji.z  + s_z*ei ;
    }
    for (unsigned int i=N_active-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        s_x  = (s_x  - p_mass[i].m * pji.x ) * ei;
        s_y  = (s_y  - p_mass[i].m * pji.y ) * ei;
        s_z  = (s_z  - p_mass[i].m * pji.z ) * ei;
        particles[i].x  = pji.x  + s_x ;
        particles[i].y  = pji.y  + s_y ;
        particles[i].z  = pji.z  + s_z ;
        eta -= p_mass[i].m;
        s_x  *= eta;
        s_y  *= eta;
        s_z  *= eta;
    }
    const double mi = 1./eta;
    particles[0].x  = s_x  * mi;
    particles[0].y  = s_y  * mi;
    particles[0].z  = s_z  * mi;
}

void reb_transformations_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active){
    double eta  = p_j[0].m;
    double s_ax  = p_j[0].ax  * eta;
    double s_ay  = p_j[0].ay  * eta;
    double s_az  = p_j[0].az  * eta;
    for (unsigned int i=N-1;i>=N_active;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        particles[i].ax  = pji.ax  + s_ax * ei;
        particles[i].ay  = pji.ay  + s_ay * ei;
        particles[i].az  = pji.az  + s_az * ei;
    }
    for (unsigned int i=N_active-1;i>0;i--){
        const struct reb_particle pji = p_j[i];
        const double ei = 1./eta;
        s_ax  = (s_ax  - p_mass[i].m * pji.ax ) * ei;
        s_ay  = (s_ay  - p_mass[i].m * pji.ay ) * ei;
        s_az  = (s_az  - p_mass[i].m * pji.az ) * ei;
        particles[i].ax  = pji.ax  + s_ax ;
        particles[i].ay  = pji.ay  + s_ay ;
        particles[i].az  = pji.az  + s_az ;
        eta -= p_mass[i].m;
        s_ax  *= eta;
        s_ay  *= eta;
        s_az  *= eta;
    }
    const double mi = 1./eta;
    particles[0].ax  = s_ax  * mi;
    particles[0].ay  = s_ay  * mi;
    particles[0].az  = s_az  * mi;
}

/******************************
 * WHDS (Hernandez)           */

void reb_transformations_inertial_to_whds_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active){
    double x0  = 0.;
    double y0  = 0.;
    double z0  = 0.;
    double vx0 = 0.;
    double vy0 = 0.;
    double vz0 = 0.;
    double m0  = 0.;
#pragma omp parallel for reduction(+:x0) reduction(+:y0) reduction(+:z0) reduction(+:vx0) reduction(+:vy0) reduction(+:vz0) reduction(+:m0)
    for (unsigned int i=0;i<N_active;i++){
        double m = particles[i].m;
        x0  += particles[i].x *m;
        y0  += particles[i].y *m;
        z0  += particles[i].z *m;
        vx0 += particles[i].vx*m;
        vy0 += particles[i].vy*m;
        vz0 += particles[i].vz*m;
        m0  += m;
    }
    p_h[0].x  = x0/m0;
    p_h[0].y  = y0/m0;
    p_h[0].z  = z0/m0;
    p_h[0].vx = vx0/m0;
    p_h[0].vy = vy0/m0;
    p_h[0].vz = vz0/m0;
    p_h[0].m = m0;
    
    m0 = particles[0].m;
#pragma omp parallel for
    for (unsigned int i=1;i<N_active;i++){
        p_h[i].x  = particles[i].x  - particles[0].x ;
        p_h[i].y  = particles[i].y  - particles[0].y ;
        p_h[i].z  = particles[i].z  - particles[0].z ;
        const double mi = particles[i].m;
        double mf = (m0+mi) / m0;
        p_h[i].vx = mf*(particles[i].vx - p_h[0].vx);
        p_h[i].vy = mf*(particles[i].vy - p_h[0].vy);
        p_h[i].vz = mf*(particles[i].vz - p_h[0].vz);
        p_h[i].m  = mi;
    }
#pragma omp parallel for
    for (unsigned int i=N_active;i<N;i++){
        p_h[i].x  = particles[i].x  - particles[0].x ;
        p_h[i].y  = particles[i].y  - particles[0].y ;
        p_h[i].z  = particles[i].z  - particles[0].z ;
        p_h[i].vx = particles[i].vx - p_h[0].vx;
        p_h[i].vy = particles[i].vy - p_h[0].vy;
        p_h[i].vz = particles[i].vz - p_h[0].vz;
        p_h[i].m  = particles[i].m;
    }
}

void reb_transformations_whds_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active){
    // Same as in heliocentric case.
    reb_transformations_democraticheliocentric_to_inertial_pos(particles, p_h, N, N_active);
}

void reb_transformations_whds_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active){
    reb_transformations_whds_to_inertial_pos(particles,p_h,N, N_active);
    const double m0 = particles[0].m;
#pragma omp parallel for
    for (unsigned int i=1;i<N_active;i++){
        const double mi = particles[i].m;
        double mf = (m0+mi) / m0;
        particles[i].vx = p_h[i].vx/mf+p_h[0].vx;
        particles[i].vy = p_h[i].vy/mf+p_h[0].vy;
        particles[i].vz = p_h[i].vz/mf+p_h[0].vz;
    }
#pragma omp parallel for
    for (unsigned int i=N_active;i<N;i++){
        particles[i].vx = p_h[i].vx+p_h[0].vx;
        particles[i].vy = p_h[i].vy+p_h[0].vy;
        particles[i].vz = p_h[i].vz+p_h[0].vz;
    }
    double vx0  = 0.;
    double vy0  = 0.;
    double vz0  = 0.;
#pragma omp parallel for reduction(+:vx0) reduction(+:vy0) reduction(+:vz0) 
    for (int i=1;i<N_active;i++){
        double m = particles[i].m;
        vx0 += p_h[i].vx*m/(m0+m);
        vy0 += p_h[i].vy*m/(m0+m);
        vz0 += p_h[i].vz*m/(m0+m);
    }
    particles[0].vx = p_h[0].vx -vx0;
    particles[0].vy = p_h[0].vy -vy0;
    particles[0].vz = p_h[0].vz -vz0;
}

/******************************
 * Democratic heliocentric.   *
 * Duncan, Levison & Lee 1998 */

void reb_transformations_inertial_to_democraticheliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active){
    double x0  = 0.;
    double y0  = 0.;
    double z0  = 0.;
    double vx0 = 0.;
    double vy0 = 0.;
    double vz0 = 0.;
    double m0  = 0.;
#pragma omp parallel for reduction(+:x0) reduction(+:y0) reduction(+:z0) reduction(+:vx0) reduction(+:vy0) reduction(+:vz0) reduction(+:m0)
    for (int i=0;i<N_active;i++){
        double m = particles[i].m;
        x0  += particles[i].x *m;
        y0  += particles[i].y *m;
        z0  += particles[i].z *m;
        vx0 += particles[i].vx*m;
        vy0 += particles[i].vy*m;
        vz0 += particles[i].vz*m;
        m0  += m;
    }
    p_h[0].x  = x0/m0;
    p_h[0].y  = y0/m0;
    p_h[0].z  = z0/m0;
    p_h[0].vx = vx0/m0;
    p_h[0].vy = vy0/m0;
    p_h[0].vz = vz0/m0;
    p_h[0].m = m0;
    
#pragma omp parallel for 
    for (unsigned int i=1;i<N;i++){
        p_h[i].x  = particles[i].x  - particles[0].x ;
        p_h[i].y  = particles[i].y  - particles[0].y ;
        p_h[i].z  = particles[i].z  - particles[0].z ;
        p_h[i].vx = particles[i].vx - p_h[0].vx;
        p_h[i].vy = particles[i].vy - p_h[0].vy;
        p_h[i].vz = particles[i].vz - p_h[0].vz;
        p_h[i].m  = particles[i].m;
    }
}

void reb_transformations_democraticheliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active){
    const double mtot = p_h[0].m;
    double x0  = 0.;
    double y0  = 0.;
    double z0  = 0.;
#pragma omp parallel for reduction(+:x0) reduction(+:y0) reduction(+:z0) 
    for (int i=1;i<N_active;i++){
        double m = p_h[i].m;
        x0 += p_h[i].x*m/mtot;
        y0 += p_h[i].y*m/mtot;
        z0 += p_h[i].z*m/mtot;
        particles[i].m = m; // in case of merger/mass change
    }
    particles[0].x  = p_h[0].x - x0;
    particles[0].y  = p_h[0].y - y0;
    particles[0].z  = p_h[0].z - z0;
#pragma omp parallel for 
    for (unsigned int i=1;i<N;i++){
        particles[i].x = p_h[i].x+particles[0].x;
        particles[i].y = p_h[i].y+particles[0].y;
        particles[i].z = p_h[i].z+particles[0].z;
    }
}

void reb_transformations_democraticheliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active){
    reb_transformations_democraticheliocentric_to_inertial_pos(particles,p_h,N,N_active);
    const double m0 = particles[0].m;
#pragma omp parallel for 
    for (unsigned int i=1;i<N;i++){
        particles[i].vx = p_h[i].vx+p_h[0].vx;
        particles[i].vy = p_h[i].vy+p_h[0].vy;
        particles[i].vz = p_h[i].vz+p_h[0].vz;
    }
    double vx0  = 0.;
    double vy0  = 0.;
    double vz0  = 0.;
#pragma omp parallel for reduction(+:vx0) reduction(+:vy0) reduction(+:vz0) 
    for (int i=1;i<N_active;i++){
        double m = particles[i].m;
        vx0 += p_h[i].vx*m/m0;
        vy0 += p_h[i].vy*m/m0;
        vz0 += p_h[i].vz*m/m0;
    }
    particles[0].vx = p_h[0].vx -vx0;
    particles[0].vy = p_h[0].vy -vy0;
    particles[0].vz = p_h[0].vz -vz0;
}

