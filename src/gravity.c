/**
 * @file     gravity.c
 * @brief     Direct gravity calculation, O(N^2).
 * @author     Hanno Rein <hanno@hanno-rein.de>
 *
 * @details     This is the crudest implementation of an N-body code
 * which sums up every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). Note that the MPI
 * implementation is not well tested and only works for very specific
 * problems. This should be resolved in the future. 
 *
 * 
 * @section LICENSE
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
#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#include "integrator_mercurius.h"
#include "integrator_trace.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

#ifdef MPI
#include "communication_mpi.h"
#endif

/**
  * @brief The function loops over all trees to call calculate_forces_for_particle_from_cell() tree to calculate forces for each particle.
  * @param r REBOUND simulation to consider
  * @param pt Index of the particle the force is calculated for.
  * @param gb Ghostbox plus position of the particle (precalculated). 
  */
static void reb_calculate_acceleration_for_particle(const struct reb_simulation* const r, const int pt, const struct reb_vec6d gb);


/**
 * Main Gravity Routine
 */
void reb_calculate_acceleration(struct reb_simulation* r){
    if (r->integrator != REB_INTEGRATOR_MERCURIUS && r->gravity == REB_GRAVITY_MERCURIUS){
        reb_simulation_warning(r,"You are using the Mercurius gravity routine with a non-Mercurius integrator. This will probably lead to unexpected behaviour. REBOUND is now setting the gravity routine back to rEB_GRAVITY_BASIC. To avoid this warning message, consider manually setting the gravity routine after changing integrators.");
        r->gravity = REB_GRAVITY_BASIC;

    }
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    const int N_active = r->N_active;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const unsigned int _gravity_ignore_terms = r->gravity_ignore_terms;
    const int _N_real   = N  - r->N_var;
    const int _N_active = ((N_active==-1)?_N_real:N_active);
    const int _testparticle_type   = r->testparticle_type;
    switch (r->gravity){
        case REB_GRAVITY_NONE: // Do nothing.
        for (int j=0; j<N; j++){
            particles[j].ax = 0;  
            particles[j].ay = 0;  
            particles[j].az = 0;  
        }  
        break;
        case REB_GRAVITY_JACOBI:
        {
            if (r->integrator != REB_INTEGRATOR_WHFAST && r->integrator != REB_INTEGRATOR_SABA ){
                reb_simulation_warning(r, "An integrator other than WHFast/SABA is being used with REB_GRAVITY_JACOBI. This is probably not correct. Use another gravity routine such as REB_GRAVITY_BASIC.");
            }
            double Rjx = 0.;
            double Rjy = 0.;
            double Rjz = 0.;
            double Mj = 0.;
            for (int j=0; j<N; j++){
                particles[j].ax = 0; 
                particles[j].ay = 0; 
                particles[j].az = 0; 
                for (int i=0; i<j+1; i++){
                    if (j>1){
                        ////////////////
                        // Jacobi Term
                        // Note: ignoring j==1 term here and below as they cancel
                        const double Qjx = particles[j].x - Rjx/Mj; 
                        const double Qjy = particles[j].y - Rjy/Mj;
                        const double Qjz = particles[j].z - Rjz/Mj;
                        const double dr = sqrt(Qjx*Qjx + Qjy*Qjy + Qjz*Qjz);
                        double dQjdri = Mj; 
                        if (i<j){
                            dQjdri = -particles[j].m; //rearranged such that m==0 does not diverge
                        }
                        const double prefact = G*dQjdri/(dr*dr*dr);
                        particles[i].ax    += prefact*Qjx;
                        particles[i].ay    += prefact*Qjy;
                        particles[i].az    += prefact*Qjz;
                    }
                    if (i!=j && (i!=0 || j!=1)){
                        ////////////////
                        // Direct Term
                        // Note: ignoring i==0 && j==1 term here and above as they cancel 
                        const double dx = particles[i].x - particles[j].x;
                        const double dy = particles[i].y - particles[j].y;
                        const double dz = particles[i].z - particles[j].z;
                        const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                        const double prefact = G /(dr*dr*dr);
                        const double prefacti = prefact*particles[i].m;
                        const double prefactj = prefact*particles[j].m;
                        
                        particles[i].ax    -= prefactj*dx;
                        particles[i].ay    -= prefactj*dy;
                        particles[i].az    -= prefactj*dz;
                        particles[j].ax    += prefacti*dx;
                        particles[j].ay    += prefacti*dy;
                        particles[j].az    += prefacti*dz;
                    }
                }
                Rjx += particles[j].m*particles[j].x;
                Rjy += particles[j].m*particles[j].y;
                Rjz += particles[j].m*particles[j].z;
                Mj += particles[j].m;
            }
        }
        break;
        case REB_GRAVITY_BASIC:
        {
            const int N_ghost_x = r->N_ghost_x;
            const int N_ghost_y = r->N_ghost_y;
            const int N_ghost_z = r->N_ghost_z;
#ifndef OPENMP // OPENMP off
            const int starti = (_gravity_ignore_terms==0)?1:2;
            const int startj = (_gravity_ignore_terms==2)?1:0;
#endif // OPENMP
#pragma omp parallel for 
            for (int i=0; i<N; i++){
                particles[i].ax = 0; 
                particles[i].ay = 0; 
                particles[i].az = 0; 
            }
            // Summing over all Ghost Boxes
            for (int gbx=-N_ghost_x; gbx<=N_ghost_x; gbx++){
            for (int gby=-N_ghost_y; gby<=N_ghost_y; gby++){
            for (int gbz=-N_ghost_z; gbz<=N_ghost_z; gbz++){
                struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                // All active particle pairs
#ifndef OPENMP // OPENMP off, do O(1/2*N^2)
                for (int i=starti; i<_N_active; i++){
                if (reb_sigint > 1) return;
                for (int j=startj; j<i; j++){
                    const double dx = (gb.x+particles[i].x) - particles[j].x;
                    const double dy = (gb.y+particles[i].y) - particles[j].y;
                    const double dz = (gb.z+particles[i].z) - particles[j].z;
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                    const double prefact = G/(_r*_r*_r);
                    const double prefactj = -prefact*particles[j].m;
                    const double prefacti = prefact*particles[i].m;
                    
                    particles[i].ax    += prefactj*dx;
                    particles[i].ay    += prefactj*dy;
                    particles[i].az    += prefactj*dz;
                    particles[j].ax    += prefacti*dx;
                    particles[j].ay    += prefacti*dy;
                    particles[j].az    += prefacti*dz;
                }
                }
#else // OPENMP on, do O(N^2)
#pragma omp parallel for
                for (int i=0; i<_N_real; i++){
                for (int j=0; j<_N_active; j++){
                    if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0) )) continue;
                    if (_gravity_ignore_terms==2 && ((j==0 || i==0) )) continue;
                    if (i==j) continue;
                    const double dx = (gb.x+particles[i].x) - particles[j].x;
                    const double dy = (gb.y+particles[i].y) - particles[j].y;
                    const double dz = (gb.z+particles[i].z) - particles[j].z;
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                    const double prefact = -G/(_r*_r*_r)*particles[j].m;
                    
                    particles[i].ax    += prefact*dx;
                    particles[i].ay    += prefact*dy;
                    particles[i].az    += prefact*dz;
                }
                }
#endif // OPENMP
                // Interactions of test particles with active particles
#ifndef OPENMP // OPENMP off
                const int startitestp = MAX(_N_active, starti);
                for (int i=startitestp; i<_N_real; i++){
                if (reb_sigint > 1) return;
                for (int j=startj; j<_N_active; j++){
                    const double dx = (gb.x+particles[i].x) - particles[j].x;
                    const double dy = (gb.y+particles[i].y) - particles[j].y;
                    const double dz = (gb.z+particles[i].z) - particles[j].z;
                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                    const double prefact = G/(_r*_r*_r);
                    const double prefactj = -prefact*particles[j].m;
                    
                    particles[i].ax    += prefactj*dx;
                    particles[i].ay    += prefactj*dy;
                    particles[i].az    += prefactj*dz;
                    if (_testparticle_type){
                        const double prefacti = prefact*particles[i].m;
                        particles[j].ax    += prefacti*dx;
                        particles[j].ay    += prefacti*dy;
                        particles[j].az    += prefacti*dz;
                    }
                }
                }
#else // OPENMP on
                if (_testparticle_type){
#pragma omp parallel for
				for (int i=0; i<_N_active; i++){
				for (int j=_N_active; j<_N_real; j++){
					if (_gravity_ignore_terms==1 && ((j==1 && i==0) )) continue;
					if (_gravity_ignore_terms==2 && ((j==0 || i==0) )) continue;
					const double dx = (gb.x+particles[i].x) - particles[j].x;
					const double dy = (gb.y+particles[i].y) - particles[j].y;
					const double dz = (gb.z+particles[i].z) - particles[j].z;
					const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
					const double prefact = -G/(_r*_r*_r)*particles[j].m;
					
					particles[i].ax    += prefact*dx;
					particles[i].ay    += prefact*dy;
					particles[i].az    += prefact*dz;
				}
				}
                }
#endif // OPENMP
            }
            }
            }
        }
        break;
        case REB_GRAVITY_COMPENSATED:
        {
            if (r->N_allocated_gravity_cs<N){
                r->gravity_cs = realloc(r->gravity_cs,N*sizeof(struct reb_vec3d));
                r->N_allocated_gravity_cs = N;
            }
            struct reb_vec3d* restrict const cs = r->gravity_cs;
#pragma omp parallel for schedule(guided)
            for (int i=0; i<_N_real; i++){
                particles[i].ax = 0.; 
                particles[i].ay = 0.; 
                particles[i].az = 0.; 
                cs[i].x = 0.;
                cs[i].y = 0.;
                cs[i].z = 0.;
            }
            // Summing over all massive particle pairs
#ifdef OPENMP
#pragma omp parallel for schedule(guided)
            for (int i=0; i<_N_active; i++){
            for (int j=0; j<_N_active; j++){
                if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                if (i==j) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefactj = -prefact*particles[j].m;
                
                {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;
                
                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
                }
            }
            }

            // Testparticles
#pragma omp parallel for schedule(guided)
            for (int i=_N_active; i<_N_real; i++){
            for (int j=0; j<_N_active; j++){
                if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefactj = -prefact*particles[j].m;
                
                {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;
                
                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
                }
            }
            }
            if (_testparticle_type){
#pragma omp parallel for schedule(guided)
                for (int j=0; j<_N_active; j++){
                for (int i=_N_active; i<_N_real; i++){
                    if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                    if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                    const double dx = particles[i].x - particles[j].x;
                    const double dy = particles[i].y - particles[j].y;
                    const double dz = particles[i].z - particles[j].z;
                    const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                    const double r = sqrt(r2);
                    const double prefact  = G/(r2*r);
                    const double prefacti = prefact*particles[i].m;
                    {
                    double ix = prefacti*dx;
                    double yx = ix - cs[j].x;
                    double tx = particles[j].ax + yx;
                    cs[j].x = (tx - particles[j].ax) - yx;
                    particles[j].ax = tx;

                    double iy = prefacti*dy;
                    double yy = iy - cs[j].y;
                    double ty = particles[j].ay + yy;
                    cs[j].y = (ty - particles[j].ay) - yy;
                    particles[j].ay = ty;
                    
                    double iz = prefacti*dz;
                    double yz = iz - cs[j].z;
                    double tz = particles[j].az + yz;
                    cs[j].z = (tz - particles[j].az) - yz;
                    particles[j].az = tz;
                    }
                }
                }
            }
#else // OPENMP
            for (int i=0; i<_N_active; i++){
            if (reb_sigint > 1) return;
            for (int j=i+1; j<_N_active; j++){
                if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefacti = prefact*particles[i].m;
                const double prefactj = -prefact*particles[j].m;
                
                {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;
                
                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
                }
                
                {
                double ix = prefacti*dx;
                double yx = ix - cs[j].x;
                double tx = particles[j].ax + yx;
                cs[j].x = (tx - particles[j].ax) - yx;
                particles[j].ax = tx;

                double iy = prefacti*dy;
                double yy = iy - cs[j].y;
                double ty = particles[j].ay + yy;
                cs[j].y = (ty - particles[j].ay) - yy;
                particles[j].ay = ty;
                
                double iz = prefacti*dz;
                double yz = iz - cs[j].z;
                double tz = particles[j].az + yz;
                cs[j].z = (tz - particles[j].az) - yz;
                particles[j].az = tz;
                }
            }
            }

            // Testparticles
            for (int i=_N_active; i<_N_real; i++){
            if (reb_sigint > 1) return;
            for (int j=0; j<_N_active; j++){
                if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double r2 = dx*dx + dy*dy + dz*dz + softening2;
                const double r = sqrt(r2);
                const double prefact  = G/(r2*r);
                const double prefactj = -prefact*particles[j].m;
                
                {
                double ix = prefactj*dx;
                double yx = ix - cs[i].x;
                double tx = particles[i].ax + yx;
                cs[i].x = (tx - particles[i].ax) - yx;
                particles[i].ax = tx;

                double iy = prefactj*dy;
                double yy = iy- cs[i].y;
                double ty = particles[i].ay + yy;
                cs[i].y = (ty - particles[i].ay) - yy;
                particles[i].ay = ty;
                
                double iz = prefactj*dz;
                double yz = iz - cs[i].z;
                double tz = particles[i].az + yz;
                cs[i].z = (tz - particles[i].az) - yz;
                particles[i].az = tz;
                }
                if (_testparticle_type){
                    const double prefacti = prefact*particles[i].m;
                    {
                    double ix = prefacti*dx;
                    double yx = ix - cs[j].x;
                    double tx = particles[j].ax + yx;
                    cs[j].x = (tx - particles[j].ax) - yx;
                    particles[j].ax = tx;

                    double iy = prefacti*dy;
                    double yy = iy - cs[j].y;
                    double ty = particles[j].ay + yy;
                    cs[j].y = (ty - particles[j].ay) - yy;
                    particles[j].ay = ty;
                    
                    double iz = prefacti*dz;
                    double yz = iz - cs[j].z;
                    double tz = particles[j].az + yz;
                    cs[j].z = (tz - particles[j].az) - yz;
                    particles[j].az = tz;
                    }
                }
            }
            }
#endif // OPENMP
        }
        break;
        case REB_GRAVITY_TREE:
        {
#pragma omp parallel for schedule(guided)
            for (int i=0; i<N; i++){
                particles[i].ax = 0; 
                particles[i].ay = 0; 
                particles[i].az = 0; 
            }
            // Summing over all Ghost Boxes
            for (int gbx=-r->N_ghost_x; gbx<=r->N_ghost_x; gbx++){
            for (int gby=-r->N_ghost_y; gby<=r->N_ghost_y; gby++){
            for (int gbz=-r->N_ghost_z; gbz<=r->N_ghost_z; gbz++){
                // Summing over all particle pairs
#pragma omp parallel for schedule(guided)
                for (int i=0; i<N; i++){
#ifndef OPENMP
                    if (reb_sigint > 1) return;
#endif // OPENMP
                    struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
                    // Precalculated shifted position
                    gb.x += particles[i].x;
                    gb.y += particles[i].y;
                    gb.z += particles[i].z;
                    reb_calculate_acceleration_for_particle(r, i, gb);
                }
            }
            }
            }
        }
        break;
        case REB_GRAVITY_MERCURIUS:
        {
            double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurius.L;
            switch (r->ri_mercurius.mode){
                case 0: // WHFAST part
                {
                    const double* const dcrit = r->ri_mercurius.dcrit;
#ifndef OPENMP
                    for (int i=0; i<_N_real; i++){
                        particles[i].ax = 0; 
                        particles[i].ay = 0; 
                        particles[i].az = 0; 
                    }
                    for (int i=2; i<_N_active; i++){
                        if (reb_sigint > 1) return;
                        for (int j=1; j<i; j++){
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[i],dcrit[j]);
                            const double L = _L(r,_r,dcritmax);
                            const double prefact = G*L/(_r*_r*_r);
                            const double prefactj = -prefact*particles[j].m;
                            const double prefacti = prefact*particles[i].m;
                            particles[i].ax    += prefactj*dx;
                            particles[i].ay    += prefactj*dy;
                            particles[i].az    += prefactj*dz;
                            particles[j].ax    += prefacti*dx;
                            particles[j].ay    += prefacti*dy;
                            particles[j].az    += prefacti*dz;
                        }
                    }
                    const int startitestp = MAX(_N_active,2);
                    for (int i=startitestp; i<_N_real; i++){
                        if (reb_sigint > 1) return;
                        for (int j=1; j<_N_active; j++){
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[i],dcrit[j]);
                            const double L = _L(r,_r,dcritmax);
                            const double prefact = G*L/(_r*_r*_r);
                            const double prefactj = -prefact*particles[j].m;
                            particles[i].ax    += prefactj*dx;
                            particles[i].ay    += prefactj*dy;
                            particles[i].az    += prefactj*dz;
                            if (_testparticle_type){
                                const double prefacti = prefact*particles[i].m;
                                particles[j].ax    += prefacti*dx;
                                particles[j].ay    += prefacti*dy;
                                particles[j].az    += prefacti*dz;
                            }
                        }
                    }
#else // OPENMP
                    particles[0].ax = 0; 
                    particles[0].ay = 0; 
                    particles[0].az = 0; 
#pragma omp parallel for schedule(guided)
                    for (int i=1; i<_N_real; i++){
                        particles[i].ax = 0; 
                        particles[i].ay = 0; 
                        particles[i].az = 0; 
                        for (int j=1; j<_N_active; j++){
                            if (i==j) continue;
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[i],dcrit[j]);
                            const double L = _L(r,_r,dcritmax);
                            const double prefact = -G*particles[j].m*L/(_r*_r*_r);
                            particles[i].ax    += prefact*dx;
                            particles[i].ay    += prefact*dy;
                            particles[i].az    += prefact*dz;
                        }
                    }
                    if (_testparticle_type){
                    for (int i=1; i<_N_active; i++){
                        for (int j=_N_active; j<_N_real; j++){
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[i],dcrit[j]);
                            const double L = _L(r,_r,dcritmax);
                            const double prefact = -G*particles[j].m*L/(_r*_r*_r);
                            particles[i].ax    += prefact*dx;
                            particles[i].ay    += prefact*dy;
                            particles[i].az    += prefact*dz;
                        }
                    }
                    }
#endif // OPENMP
                }
                break;
                case 1: // IAS15 part
                {
                    const double m0 = r->particles[0].m;
                    const double* const dcrit = r->ri_mercurius.dcrit;
                    const int encounter_N = r->ri_mercurius.encounter_N;
                    const int encounter_N_active = r->ri_mercurius.encounter_N_active;
                    int* map = r->ri_mercurius.encounter_map;
#ifndef OPENMP
                    particles[0].ax = 0; // map[0] is always 0 
                    particles[0].ay = 0; 
                    particles[0].az = 0; 
                    // Acceleration due to star
                    for (int i=1; i<encounter_N; i++){
                        int mi = map[i];
                        const double x = particles[mi].x;
                        const double y = particles[mi].y;
                        const double z = particles[mi].z;
                        const double _r = sqrt(x*x + y*y + z*z + softening2);
                        double prefact = -G/(_r*_r*_r)*m0;
                        particles[mi].ax    = prefact*x;
                        particles[mi].ay    = prefact*y;
                        particles[mi].az    = prefact*z;
                    }
                    // We're in a heliocentric coordinate system.
                    // The star feels no acceleration
                    // Interactions between active-active
                    for (int i=2; i<encounter_N_active; i++){
                        int mi = map[i];
                        for (int j=1; j<i; j++){
                            int mj = map[j];
                            const double dx = particles[mi].x - particles[mj].x;
                            const double dy = particles[mi].y - particles[mj].y;
                            const double dz = particles[mi].z - particles[mj].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                            const double L = _L(r,_r,dcritmax);
                            double prefact = G*(1.-L)/(_r*_r*_r);
                            double prefactj = -prefact*particles[mj].m;
                            double prefacti = prefact*particles[mi].m;
                            particles[mi].ax    += prefactj*dx;
                            particles[mi].ay    += prefactj*dy;
                            particles[mi].az    += prefactj*dz;
                            particles[mj].ax    += prefacti*dx;
                            particles[mj].ay    += prefacti*dy;
                            particles[mj].az    += prefacti*dz;
                        }
                    }
                    // Interactions between active-testparticle
                    const int startitestp = MAX(encounter_N_active,2);
                    for (int i=startitestp; i<encounter_N; i++){
                        int mi = map[i];
                        for (int j=1; j<encounter_N_active; j++){
                            int mj = map[j];
                            const double dx = particles[mi].x - particles[mj].x;
                            const double dy = particles[mi].y - particles[mj].y;
                            const double dz = particles[mi].z - particles[mj].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                            const double L = _L(r,_r,dcritmax);
                            double prefact = G*(1.-L)/(_r*_r*_r);
                            double prefactj = -prefact*particles[mj].m;
                            particles[mi].ax    += prefactj*dx;
                            particles[mi].ay    += prefactj*dy;
                            particles[mi].az    += prefactj*dz;
                            if (_testparticle_type){
                                double prefacti = prefact*particles[mi].m;
                                particles[mj].ax    += prefacti*dx;
                                particles[mj].ay    += prefacti*dy;
                                particles[mj].az    += prefacti*dz;
                            }
                        }
                    }
#else // OPENMP
                    particles[0].ax = 0; // map[0] is always 0 
                    particles[0].ay = 0; 
                    particles[0].az = 0; 
                    // We're in a heliocentric coordinate system.
                    // The star feels no acceleration
#pragma omp parallel for schedule(guided)
                    for (int i=1; i<encounter_N; i++){
                        int mi = map[i];
                        particles[mi].ax = 0; 
                        particles[mi].ay = 0; 
                        particles[mi].az = 0; 
                        // Acceleration due to star
                        const double x = particles[mi].x;
                        const double y = particles[mi].y;
                        const double z = particles[mi].z;
                        const double _r = sqrt(x*x + y*y + z*z + softening2);
                        double prefact = -G/(_r*_r*_r)*m0;
                        particles[mi].ax    += prefact*x;
                        particles[mi].ay    += prefact*y;
                        particles[mi].az    += prefact*z;
                        for (int j=1; j<encounter_N_active; j++){
                            if (i==j) continue;
                            int mj = map[j];
                            const double dx = x - particles[mj].x;
                            const double dy = y - particles[mj].y;
                            const double dz = z - particles[mj].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                            const double L = _L(r,_r,dcritmax);
                            double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
                            particles[mi].ax    += prefact*dx;
                            particles[mi].ay    += prefact*dy;
                            particles[mi].az    += prefact*dz;
                        }
                    }
                    if (_testparticle_type){
#pragma omp parallel for schedule(guided)
                    for (int i=1; i<encounter_N_active; i++){
                        int mi = map[i];
                        const double x = particles[mi].x;
                        const double y = particles[mi].y;
                        const double z = particles[mi].z;
                        for (int j=encounter_N_active; j<encounter_N; j++){
                            int mj = map[j];
                            const double dx = x - particles[mj].x;
                            const double dy = y - particles[mj].y;
                            const double dz = z - particles[mj].z;
                            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                            const double L = _L(r,_r,dcritmax);
                            double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
                            particles[mi].ax    += prefact*dx;
                            particles[mi].ay    += prefact*dy;
                            particles[mi].az    += prefact*dz;
                        }
                    }
                    }
#endif // OPENMP
                }
                break;
                case 2: // Skip WHFAST part because of synchronization
                break;
              }
            }
              break;
                case REB_GRAVITY_TRACE:
                {
                    switch (r->ri_trace.mode){
                        case REB_TRACE_MODE_INTERACTION: // Interaction step
                        {
        #ifndef OPENMP
                            for (int i=0; i<_N_real; i++){
                                particles[i].ax = 0;
                                particles[i].ay = 0;
                                particles[i].az = 0;
                            }
                            for (int i=2; i<_N_active; i++){
                                if (reb_sigint > 1) return;
                                for (int j=1; j<i; j++){
                                    if (r->ri_trace.current_Ks[j*N+i]) continue;
                                    const double dx = particles[i].x - particles[j].x;
                                    const double dy = particles[i].y - particles[j].y;
                                    const double dz = particles[i].z - particles[j].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    const double prefact = G / (_r*_r*_r);
                                    const double prefactj = -prefact*particles[j].m;
                                    const double prefacti = prefact*particles[i].m;
                                    particles[i].ax    += prefactj*dx;
                                    particles[i].ay    += prefactj*dy;
                                    particles[i].az    += prefactj*dz;
                                    particles[j].ax    += prefacti*dx;
                                    particles[j].ay    += prefacti*dy;
                                    particles[j].az    += prefacti*dz;
                                }
                            }
                            const int startitestp = MAX(_N_active,2);
                            for (int i=startitestp; i<_N_real; i++){
                                if (reb_sigint > 1) return;
                                for (int j=1; j<_N_active; j++){
                                    if (r->ri_trace.current_Ks[j*N+i]) continue;
                                    const double dx = particles[i].x - particles[j].x;
                                    const double dy = particles[i].y - particles[j].y;
                                    const double dz = particles[i].z - particles[j].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    const double prefact = G / (_r*_r*_r);
                                    const double prefactj = -prefact*particles[j].m;
                                    particles[i].ax    += prefactj*dx;
                                    particles[i].ay    += prefactj*dy;
                                    particles[i].az    += prefactj*dz;
                                    if (_testparticle_type){
                                        const double prefacti = prefact*particles[i].m;
                                        particles[j].ax    += prefacti*dx;
                                        particles[j].ay    += prefacti*dy;
                                        particles[j].az    += prefacti*dz;
                                    }
                                }
                            }

        #else // OPENMP
                            particles[0].ax = 0;
                            particles[0].ay = 0;
                            particles[0].az = 0;
        #pragma omp parallel for schedule(guided)
                            for (int i=1; i<_N_real; i++){
                                particles[i].ax = 0;
                                particles[i].ay = 0;
                                particles[i].az = 0;
                                for (int j=1; j<_N_active; j++){
                                    if (i==j) continue;
                                    if (r->ri_trace.current_Ks[j*N+i]) continue;
                                    const double dx = particles[i].x - particles[j].x;
                                    const double dy = particles[i].y - particles[j].y;
                                    const double dz = particles[i].z - particles[j].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    const double prefact = -G*particles[j].m/(_r*_r*_r);
                                    particles[i].ax    += prefact*dx;
                                    particles[i].ay    += prefact*dy;
                                    particles[i].az    += prefact*dz;
                                }
                            }
                            if (_testparticle_type){
                            for (int i=1; i<_N_active; i++){
                                for (int j=_N_active; j<_N_real; j++){
                                    if (r->ri_trace.current_Ks[j*N+i]) continue;
                                    const double dx = particles[i].x - particles[j].x;
                                    const double dy = particles[i].y - particles[j].y;
                                    const double dz = particles[i].z - particles[j].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    const double prefact = -G*particles[j].m/(_r*_r*_r);
                                    particles[i].ax    += prefact*dx;
                                    particles[i].ay    += prefact*dy;
                                    particles[i].az    += prefact*dz;
                                }
                            }
                            }
        #endif // OPENMP
                        }
                        break;
                        case REB_TRACE_MODE_KEPLER: // BS part
                        // Kepler Step
                        {
                            const double m0 = r->particles[0].m;
                            const int encounter_N = r->ri_trace.encounter_N;
                            const int encounter_N_active = r->ri_trace.encounter_N_active;
                            int* map = r->ri_trace.encounter_map;
        #ifndef OPENMP
                            particles[0].ax = 0; // map[0] is always 0
                            particles[0].ay = 0;
                            particles[0].az = 0;

                            // Acceleration due to star
                            for (int i=1; i<encounter_N; i++){
                                int mi = map[i];
                                const double x = particles[mi].x;
                                const double y = particles[mi].y;
                                const double z = particles[mi].z;
                                const double _r = sqrt(x*x + y*y + z*z + softening2);
                                double prefact = -G * m0 / (_r*_r*_r);
                                particles[mi].ax    = prefact*x;
                                particles[mi].ay    = prefact*y;
                                particles[mi].az    = prefact*z;
                            }

                            // We're in a heliocentric coordinate system.
                            // The star feels no acceleration
                            // Interactions between active-active
                            if (encounter_N_active > 2){ // if two or less, no active-active planets
                                for (int i=2; i<encounter_N_active; i++){
                                    int mi = map[i];
                                    for (int j=1; j<i; j++){
                                        int mj = map[j];
                                        if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                                        const double dx = particles[mi].x - particles[mj].x;
                                        const double dy = particles[mi].y - particles[mj].y;
                                        const double dz = particles[mi].z - particles[mj].z;
                                        const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                        double prefact = G/(_r*_r*_r);
                                        double prefactj = -prefact*particles[mj].m;
                                        double prefacti = prefact*particles[mi].m;

                                        particles[mi].ax    += prefactj*dx;
                                        particles[mi].ay    += prefactj*dy;
                                        particles[mi].az    += prefactj*dz;
                                        particles[mj].ax    += prefacti*dx;
                                        particles[mj].ay    += prefacti*dy;
                                        particles[mj].az    += prefacti*dz;
                                    }
                                }
                            }

                            // Interactions between active-testparticle
                            const int startitestp = MAX(encounter_N_active,2);
                            for (int i=startitestp; i<encounter_N; i++){
                                int mi = map[i];
                                for (int j=1; j<encounter_N_active; j++){
                                    int mj = map[j];
                                    if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                                    const double dx = particles[mi].x - particles[mj].x;
                                    const double dy = particles[mi].y - particles[mj].y;
                                    const double dz = particles[mi].z - particles[mj].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    double prefact = G/(_r*_r*_r);
                                    double prefactj = -prefact*particles[mj].m;
                                    particles[mi].ax    += prefactj*dx;
                                    particles[mi].ay    += prefactj*dy;
                                    particles[mi].az    += prefactj*dz;

                                    if (_testparticle_type){
                                        double prefacti = prefact*particles[mi].m;
                                        particles[mj].ax    += prefacti*dx;
                                        particles[mj].ay    += prefacti*dy;
                                        particles[mj].az    += prefacti*dz;
                                    }
                                }
                            }
        #else // OPENMP
                            particles[0].ax = 0; // map[0] is always 0
                            particles[0].ay = 0;
                            particles[0].az = 0;
                            // We're in a heliocentric coordinate system.
                            // The star feels no acceleration
#pragma omp parallel for schedule(guided)
                            for (int i=1; i<encounter_N; i++){
                                int mi = map[i];
                                particles[mi].ax = 0;
                                particles[mi].ay = 0;
                                particles[mi].az = 0;
                                // Acceleration due to star
                                const double x = particles[mi].x;
                                const double y = particles[mi].y;
                                const double z = particles[mi].z;
                                const double _r = sqrt(x*x + y*y + z*z + softening2);
                                double prefact = -G/(_r*_r*_r)*m0;
                                particles[mi].ax    += prefact*x;
                                particles[mi].ay    += prefact*y;
                                particles[mi].az    += prefact*z;
                                for (int j=1; j<encounter_N_active; j++){
                                    if (i==j) continue;
                                    int mj = map[j];
                                    if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                                    const double dx = x - particles[mj].x;
                                    const double dy = y - particles[mj].y;
                                    const double dz = z - particles[mj].z;
                                    const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                    double prefact = -G*particles[mj].m/(_r*_r*_r);
                                    particles[mi].ax    += prefact*dx;
                                    particles[mi].ay    += prefact*dy;
                                    particles[mi].az    += prefact*dz;
                                }
                            }
                            if (_testparticle_type){
#pragma omp parallel for schedule(guided)
                                for (int i=1; i<encounter_N_active; i++){
                                    int mi = map[i];
                                    const double x = particles[mi].x;
                                    const double y = particles[mi].y;
                                    const double z = particles[mi].z;
                                    for (int j=encounter_N_active; j<encounter_N; j++){
                                        int mj = map[j];
                                        if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                                        const double dx = x - particles[mj].x;
                                        const double dy = y - particles[mj].y;
                                        const double dz = z - particles[mj].z;
                                        const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                                        double prefact = -G*particles[mj].m/(_r*_r*_r);
                                        particles[mi].ax    += prefact*dx;
                                        particles[mi].ay    += prefact*dy;
                                        particles[mi].az    += prefact*dz;
                                    }
                                }
                            }
        #endif // OPENMP
                        }
                        break;
                        case REB_TRACE_MODE_NONE: // In-between steps. Do not calculate anything. 
                            break;
                        default:
                            reb_simulation_error(r, "TRACE mode not supported in gravity.c");
                            break;
                    }
        }
        break;
        default:
            reb_exit("Gravity calculation not yet implemented.");
    }

}

void reb_calculate_acceleration_var(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    const double G = r->G;
    const unsigned int _gravity_ignore_terms = r->gravity_ignore_terms;
    const int _testparticle_type   = r->testparticle_type;
    const int N = r->N;
    const int _N_real   = N - r->N_var;
    const int _N_active = ((r->N_active==-1)?_N_real:r->N_active);
    const int starti = (r->gravity_ignore_terms==0)?1:2;
    const int startj = (r->gravity_ignore_terms==2)?1:0;
    switch (r->gravity){
        case REB_GRAVITY_NONE: // Do nothing.
        break;
        case REB_GRAVITY_COMPENSATED:
        {
            struct reb_vec3d* restrict const cs = r->gravity_cs;
#pragma omp parallel for schedule(guided)
            for (int i=_N_real; i<N; i++){
                cs[i].x = 0.;
                cs[i].y = 0.;
                cs[i].z = 0.;
            }
        }
        // Fallthrough is on purpose.
        case REB_GRAVITY_BASIC:
            for (int v=0;v<r->N_var_config;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                if (vc.order==1){
                    //////////////////
                    /// 1st order  ///
                    //////////////////
                    struct reb_particle* const particles_var1 = particles + vc.index;
                    if (vc.testparticle<0){
                        for (int i=0; i<_N_real; i++){
                            particles_var1[i].ax = 0.; 
                            particles_var1[i].ay = 0.; 
                            particles_var1[i].az = 0.; 
                        }
                        for (int i=starti; i<_N_active; i++){
                        for (int j=startj; j<i; j++){
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double r2 = dx*dx + dy*dy + dz*dz;
                            const double _r  = sqrt(r2);
                            const double r3inv = 1./(r2*_r);
                            const double r5inv = 3.*r3inv/r2;
                            const double ddx = particles_var1[i].x - particles_var1[j].x;
                            const double ddy = particles_var1[i].y - particles_var1[j].y;
                            const double ddz = particles_var1[i].z - particles_var1[j].z;
                            const double Gmi = G * particles[i].m;
                            const double Gmj = G * particles[j].m;

                            // Variational equations
                            const double dxdx = dx*dx*r5inv - r3inv;
                            const double dydy = dy*dy*r5inv - r3inv;
                            const double dzdz = dz*dz*r5inv - r3inv;
                            const double dxdy = dx*dy*r5inv;
                            const double dxdz = dx*dz*r5inv;
                            const double dydz = dy*dz*r5inv;
                            const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                            const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                            const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                            // Variational mass contributions
                            const double dGmi = G*particles_var1[i].m;
                            const double dGmj = G*particles_var1[j].m;

                            particles_var1[i].ax += Gmj * dax - dGmj*r3inv*dx;
                            particles_var1[i].ay += Gmj * day - dGmj*r3inv*dy;
                            particles_var1[i].az += Gmj * daz - dGmj*r3inv*dz;

                            particles_var1[j].ax -= Gmi * dax - dGmi*r3inv*dx;
                            particles_var1[j].ay -= Gmi * day - dGmi*r3inv*dy;
                            particles_var1[j].az -= Gmi * daz - dGmi*r3inv*dz; 
                        }
                        }
                        for (int i=_N_active; i<_N_real; i++){
                        for (int j=startj; j<_N_active; j++){
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double r2 = dx*dx + dy*dy + dz*dz;
                            const double _r  = sqrt(r2);
                            const double r3inv = 1./(r2*_r);
                            const double r5inv = 3.*r3inv/r2;
                            const double ddx = particles_var1[i].x - particles_var1[j].x;
                            const double ddy = particles_var1[i].y - particles_var1[j].y;
                            const double ddz = particles_var1[i].z - particles_var1[j].z;
                            const double Gmi = G * particles[i].m;
                            const double Gmj = G * particles[j].m;

                            // Variational equations
                            const double dxdx = dx*dx*r5inv - r3inv;
                            const double dydy = dy*dy*r5inv - r3inv;
                            const double dzdz = dz*dz*r5inv - r3inv;
                            const double dxdy = dx*dy*r5inv;
                            const double dxdz = dx*dz*r5inv;
                            const double dydz = dy*dz*r5inv;
                            const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                            const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                            const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                            // Variational mass contributions
                            const double dGmi = G*particles_var1[i].m;
                            const double dGmj = G*particles_var1[j].m;

                            particles_var1[i].ax += Gmj * dax - dGmj*r3inv*dx;
                            particles_var1[i].ay += Gmj * day - dGmj*r3inv*dy;
                            particles_var1[i].az += Gmj * daz - dGmj*r3inv*dz;
                            if (_testparticle_type){
                                // Warning! This does not make sense when the mass is varied!
                                particles_var1[j].ax -= Gmi * dax - dGmi*r3inv*dx;
                                particles_var1[j].ay -= Gmi * day - dGmi*r3inv*dy;
                                particles_var1[j].az -= Gmi * daz - dGmi*r3inv*dz; 
                            }
                        }
                        }
                    }else{ //testparticle
                        int i = vc.testparticle;
                        particles_var1[0].ax = 0.; 
                        particles_var1[0].ay = 0.; 
                        particles_var1[0].az = 0.; 
                        for (int j=0; j<_N_real; j++){
                            if (i==j) continue;
                            if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                            if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double r2 = dx*dx + dy*dy + dz*dz;
                            const double _r  = sqrt(r2);
                            const double r3inv = 1./(r2*_r);
                            const double r5inv = 3.*r3inv/r2;
                            const double ddx = particles_var1[0].x;
                            const double ddy = particles_var1[0].y;
                            const double ddz = particles_var1[0].z;
                            const double Gmj = G * particles[j].m;

                            // Variational equations
                            const double dxdx = dx*dx*r5inv - r3inv;
                            const double dydy = dy*dy*r5inv - r3inv;
                            const double dzdz = dz*dz*r5inv - r3inv;
                            const double dxdy = dx*dy*r5inv;
                            const double dxdz = dx*dz*r5inv;
                            const double dydz = dy*dz*r5inv;
                            const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
                            const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
                            const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

                            // No variational mass contributions for test particles!

                            particles_var1[0].ax += Gmj * dax;
                            particles_var1[0].ay += Gmj * day;
                            particles_var1[0].az += Gmj * daz;

                        }
                    }
                }else if (vc.order==2){
                    if (_testparticle_type){
                        reb_simulation_error(r,"testparticletype=1 not implemented for second order variational equations.");
                    }
                    //////////////////
                    /// 2nd order  ///
                    //////////////////
                    struct reb_particle* const particles_var2 = particles + vc.index;
                    struct reb_particle* const particles_var1a = particles + vc.index_1st_order_a;
                    struct reb_particle* const particles_var1b = particles + vc.index_1st_order_b;
                    if (vc.testparticle<0){
                        for (int i=0; i<_N_real; i++){
                            particles_var2[i].ax = 0.; 
                            particles_var2[i].ay = 0.; 
                            particles_var2[i].az = 0.; 
                        }
                        for (int i=0; i<_N_real; i++){
                        for (int j=i+1; j<_N_real; j++){
                            // TODO: Need to implement WH skipping
                            //if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                            //if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double r2 = dx*dx + dy*dy + dz*dz;
                            const double r  = sqrt(r2);
                            const double r3inv = 1./(r2*r);
                            const double r5inv = r3inv/r2;
                            const double r7inv = r5inv/r2;
                            const double ddx = particles_var2[i].x - particles_var2[j].x;
                            const double ddy = particles_var2[i].y - particles_var2[j].y;
                            const double ddz = particles_var2[i].z - particles_var2[j].z;
                            const double Gmi = G * particles[i].m;
                            const double Gmj = G * particles[j].m;
                            const double ddGmi = G*particles_var2[i].m;
                            const double ddGmj = G*particles_var2[j].m;
                            
                            // Variational equations
                            // delta^(2) terms
                            double dax =         ddx * ( 3.*dx*dx*r5inv - r3inv )
                                       + ddy * ( 3.*dx*dy*r5inv )
                                       + ddz * ( 3.*dx*dz*r5inv );
                            double day =         ddx * ( 3.*dy*dx*r5inv )
                                       + ddy * ( 3.*dy*dy*r5inv - r3inv )
                                       + ddz * ( 3.*dy*dz*r5inv );
                            double daz =         ddx * ( 3.*dz*dx*r5inv )
                                       + ddy * ( 3.*dz*dy*r5inv )
                                       + ddz * ( 3.*dz*dz*r5inv - r3inv );
                            
                            // delta^(1) delta^(1) terms
                            const double dk1dx = particles_var1a[i].x - particles_var1a[j].x;
                            const double dk1dy = particles_var1a[i].y - particles_var1a[j].y;
                            const double dk1dz = particles_var1a[i].z - particles_var1a[j].z;
                            const double dk2dx = particles_var1b[i].x - particles_var1b[j].x;
                            const double dk2dy = particles_var1b[i].y - particles_var1b[j].y;
                            const double dk2dz = particles_var1b[i].z - particles_var1b[j].z;

                            const double rdk1 =  dx*dk1dx + dy*dk1dy + dz*dk1dz;
                            const double rdk2 =  dx*dk2dx + dy*dk2dy + dz*dk2dz;
                            const double dk1dk2 =  dk1dx*dk2dx + dk1dy*dk2dy + dk1dz*dk2dz;
                            dax     +=        3.* r5inv * dk2dx * rdk1
                                    + 3.* r5inv * dk1dx * rdk2
                                    + 3.* r5inv    * dx * dk1dk2  
                                        - 15.      * dx * r7inv * rdk1 * rdk2;
                            day     +=        3.* r5inv * dk2dy * rdk1
                                    + 3.* r5inv * dk1dy * rdk2
                                    + 3.* r5inv    * dy * dk1dk2  
                                        - 15.      * dy * r7inv * rdk1 * rdk2;
                            daz     +=        3.* r5inv * dk2dz * rdk1
                                    + 3.* r5inv * dk1dz * rdk2
                                    + 3.* r5inv    * dz * dk1dk2  
                                        - 15.      * dz * r7inv * rdk1 * rdk2;
                            
                            const double dk1Gmi = G * particles_var1a[i].m;
                            const double dk1Gmj = G * particles_var1a[j].m;
                            const double dk2Gmi = G * particles_var1b[i].m;
                            const double dk2Gmj = G * particles_var1b[j].m;

                            particles_var2[i].ax += Gmj * dax 
                                - ddGmj*r3inv*dx 
                                - dk2Gmj*r3inv*dk1dx + 3.*dk2Gmj*r5inv*dx*rdk1
                                - dk1Gmj*r3inv*dk2dx + 3.*dk1Gmj*r5inv*dx*rdk2;
                            particles_var2[i].ay += Gmj * day 
                                - ddGmj*r3inv*dy
                                - dk2Gmj*r3inv*dk1dy + 3.*dk2Gmj*r5inv*dy*rdk1
                                - dk1Gmj*r3inv*dk2dy + 3.*dk1Gmj*r5inv*dy*rdk2;
                            particles_var2[i].az += Gmj * daz 
                                - ddGmj*r3inv*dz
                                - dk2Gmj*r3inv*dk1dz + 3.*dk2Gmj*r5inv*dz*rdk1
                                - dk1Gmj*r3inv*dk2dz + 3.*dk1Gmj*r5inv*dz*rdk2;
                                                                                 
                            particles_var2[j].ax -= Gmi * dax 
                                - ddGmi*r3inv*dx
                                - dk2Gmi*r3inv*dk1dx + 3.*dk2Gmi*r5inv*dx*rdk1
                                - dk1Gmi*r3inv*dk2dx + 3.*dk1Gmi*r5inv*dx*rdk2;
                            particles_var2[j].ay -= Gmi * day 
                                - ddGmi*r3inv*dy
                                - dk2Gmi*r3inv*dk1dy + 3.*dk2Gmi*r5inv*dy*rdk1
                                - dk1Gmi*r3inv*dk2dy + 3.*dk1Gmi*r5inv*dy*rdk2;
                            particles_var2[j].az -= Gmi * daz 
                                - ddGmi*r3inv*dz
                                - dk2Gmi*r3inv*dk1dz + 3.*dk2Gmi*r5inv*dz*rdk1
                                - dk1Gmi*r3inv*dk2dz + 3.*dk1Gmi*r5inv*dz*rdk2;
                        }
                        }
                    }else{ //testparticle
                        int i = vc.testparticle;
                        particles_var2[0].ax = 0.; 
                        particles_var2[0].ay = 0.; 
                        particles_var2[0].az = 0.; 
                        for (int j=0; j<_N_real; j++){
                            if (i==j) continue;
                            // TODO: Need to implement WH skipping
                            //if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                            //if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
                            const double dx = particles[i].x - particles[j].x;
                            const double dy = particles[i].y - particles[j].y;
                            const double dz = particles[i].z - particles[j].z;
                            const double r2 = dx*dx + dy*dy + dz*dz;
                            const double r  = sqrt(r2);
                            const double r3inv = 1./(r2*r);
                            const double r5inv = r3inv/r2;
                            const double r7inv = r5inv/r2;
                            const double ddx = particles_var2[0].x;
                            const double ddy = particles_var2[0].y;
                            const double ddz = particles_var2[0].z;
                            const double Gmj = G * particles[j].m;
                            
                            // Variational equations
                            // delta^(2) terms
                            double dax =         ddx * ( 3.*dx*dx*r5inv - r3inv )
                                       + ddy * ( 3.*dx*dy*r5inv )
                                       + ddz * ( 3.*dx*dz*r5inv );
                            double day =         ddx * ( 3.*dy*dx*r5inv )
                                       + ddy * ( 3.*dy*dy*r5inv - r3inv )
                                       + ddz * ( 3.*dy*dz*r5inv );
                            double daz =         ddx * ( 3.*dz*dx*r5inv )
                                       + ddy * ( 3.*dz*dy*r5inv )
                                       + ddz * ( 3.*dz*dz*r5inv - r3inv );
                            
                            // delta^(1) delta^(1) terms
                            const double dk1dx = particles_var1a[0].x;
                            const double dk1dy = particles_var1a[0].y;
                            const double dk1dz = particles_var1a[0].z;
                            const double dk2dx = particles_var1b[0].x;
                            const double dk2dy = particles_var1b[0].y;
                            const double dk2dz = particles_var1b[0].z;

                            const double rdk1 =  dx*dk1dx + dy*dk1dy + dz*dk1dz;
                            const double rdk2 =  dx*dk2dx + dy*dk2dy + dz*dk2dz;
                            const double dk1dk2 =  dk1dx*dk2dx + dk1dy*dk2dy + dk1dz*dk2dz;
                            dax     +=        3.* r5inv * dk2dx * rdk1
                                    + 3.* r5inv * dk1dx * rdk2
                                    + 3.* r5inv    * dx * dk1dk2  
                                        - 15.      * dx * r7inv * rdk1 * rdk2;
                            day     +=        3.* r5inv * dk2dy * rdk1
                                    + 3.* r5inv * dk1dy * rdk2
                                    + 3.* r5inv    * dy * dk1dk2  
                                        - 15.      * dy * r7inv * rdk1 * rdk2;
                            daz     +=        3.* r5inv * dk2dz * rdk1
                                    + 3.* r5inv * dk1dz * rdk2
                                    + 3.* r5inv    * dz * dk1dk2  
                                        - 15.      * dz * r7inv * rdk1 * rdk2;
                            
                            // No variational mass contributions for test particles!

                            particles_var2[0].ax += Gmj * dax; 
                            particles_var2[0].ay += Gmj * day;
                            particles_var2[0].az += Gmj * daz;
                        }
                    }
                }
            }
            break;
        default:
            reb_exit("Variational gravity calculation not yet implemented.");
    }

}

void reb_calculate_and_apply_jerk(struct reb_simulation* r, const double v){
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    const int N_active = r->N_active;
    const double G = r->G;
    const int _N_real   = N  - r->N_var;
    const int _N_active = ((N_active==-1)?_N_real:N_active);
    const int _testparticle_type   = r->testparticle_type;
    const int starti = (r->gravity_ignore_terms==0)?1:2;
    const int startj = (r->gravity_ignore_terms==2)?1:0;
    switch (r->gravity){
        case REB_GRAVITY_NONE: // Do nothing.
        break;
        case REB_GRAVITY_BASIC:
            // All interactions between active particles
#pragma omp parallel for
            for (int i=starti; i<_N_active; i++){
#ifndef OPENMP
                if (reb_sigint > 1) return;
#endif // OPENMP
                for (int j=startj; j<i; j++){
                    const double dx = particles[i].x - particles[j].x; 
                    const double dy = particles[i].y - particles[j].y; 
                    const double dz = particles[i].z - particles[j].z; 
                    
                    const double dax = particles[i].ax - particles[j].ax; 
                    const double day = particles[i].ay - particles[j].ay; 
                    const double daz = particles[i].az - particles[j].az; 

                    const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                    const double alphasum = dax*dx+day*dy+daz*dz;
                    const double prefact2 = 2.*v*G /(dr*dr*dr);
                    const double prefact2i = prefact2*particles[j].m;
                    const double prefact2j = prefact2*particles[i].m;
                    const double prefact1 = alphasum*prefact2/dr *3./dr;
                    const double prefact1i = prefact1*particles[j].m;
                    const double prefact1j = prefact1*particles[i].m;
                    particles[i].vx    += dx*prefact1i - dax*prefact2i;
                    particles[i].vy    += dy*prefact1i - day*prefact2i;
                    particles[i].vz    += dz*prefact1i - daz*prefact2i;
                    particles[j].vx    += dax*prefact2j - dx*prefact1j;
                    particles[j].vy    += day*prefact2j - dy*prefact1j;
                    particles[j].vz    += daz*prefact2j - dz*prefact1j;
                }
            }
            // Interactions between active particles and test particles
#pragma omp parallel for
            for (int i=_N_active; i<_N_real; i++){
#ifndef OPENMP
                if (reb_sigint > 1) return;
#endif // OPENMP
                for (int j=startj; j<i; j++){
                    const double dx = particles[i].x - particles[j].x; 
                    const double dy = particles[i].y - particles[j].y; 
                    const double dz = particles[i].z - particles[j].z; 
                    
                    const double dax = particles[i].ax - particles[j].ax; 
                    const double day = particles[i].ay - particles[j].ay; 
                    const double daz = particles[i].az - particles[j].az; 

                    const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                    const double alphasum = dax*dx+day*dy+daz*dz;
                    const double prefact2 = 2.*v*G /(dr*dr*dr);
                    const double prefact1 = alphasum*prefact2/dr *3./dr;
                    const double prefact1i = prefact1*particles[j].m;
                    const double prefact2i = prefact2*particles[j].m;
                    particles[i].vx    += dx*prefact1i - dax*prefact2i;
                    particles[i].vy    += dy*prefact1i - day*prefact2i;
                    particles[i].vz    += dz*prefact1i - daz*prefact2i;
                    if (_testparticle_type){
                        const double prefact1j = prefact1*particles[i].m;
                        const double prefact2j = prefact2*particles[i].m;
                        particles[j].vx    += dax*prefact2j - dx*prefact1j;
                        particles[j].vy    += day*prefact2j - dy*prefact1j;
                        particles[j].vz    += daz*prefact2j - dz*prefact1j;
                    }
                }
            }
            break;
        default:
            reb_simulation_error(r,"Jerk calculation only supported for BASIC gravity routine.");
        break;
    }
}

// Helper routines for REB_GRAVITY_TREE


/**
  * @brief The function calls itself recursively using cell breaking criterion to check whether it can use center of mass (and mass quadrupole tensor) to calculate forces.
  * Calculate the acceleration for a particle from a given cell and all its daughter cells.
  *
  * @param r REBOUND simulation to consider
  * @param pt Index of the particle the force is calculated for.
  * @param node Pointer to the cell the force is calculated from.
  * @param gb Ghostbox plus position of the particle (precalculated). 
  */
static void reb_calculate_acceleration_for_particle_from_cell(const struct reb_simulation* const r, const int pt, const struct reb_treecell *node, const struct reb_vec6d gb);

static void reb_calculate_acceleration_for_particle(const struct reb_simulation* const r, const int pt, const struct reb_vec6d gb) {
    for(int i=0;i<r->N_root;i++){
        struct reb_treecell* node = r->tree_root[i];
        if (node!=NULL){
            reb_calculate_acceleration_for_particle_from_cell(r, pt, node, gb);
        }
    }
}

static void reb_calculate_acceleration_for_particle_from_cell(const struct reb_simulation* r, const int pt, const struct reb_treecell *node, const struct reb_vec6d gb) {
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    struct reb_particle* const particles = r->particles;
    const double dx = gb.x - node->mx;
    const double dy = gb.y - node->my;
    const double dz = gb.z - node->mz;
    const double r2 = dx*dx + dy*dy + dz*dz;
    if ( node->pt < 0 ) { // Not a leaf
        if ( node->w*node->w > r->opening_angle2*r2 ){
            for (int o=0; o<8; o++) {
                if (node->oct[o] != NULL) {
                    reb_calculate_acceleration_for_particle_from_cell(r, pt, node->oct[o], gb);
                }
            }
        } else {
            double _r = sqrt(r2 + softening2);
            double prefact = -G/(_r*_r*_r)*node->m;
#ifdef QUADRUPOLE
            double qprefact = G/(_r*_r*_r*_r*_r);
            particles[pt].ax += qprefact*(dx*node->mxx + dy*node->mxy + dz*node->mxz); 
            particles[pt].ay += qprefact*(dx*node->mxy + dy*node->myy + dz*node->myz); 
            particles[pt].az += qprefact*(dx*node->mxz + dy*node->myz + dz*node->mzz); 
            double mrr     = dx*dx*node->mxx     + dy*dy*node->myy     + dz*dz*node->mzz
                    + 2.*dx*dy*node->mxy     + 2.*dx*dz*node->mxz     + 2.*dy*dz*node->myz; 
            qprefact *= -5.0/(2.0*_r*_r)*mrr;
            particles[pt].ax += (qprefact + prefact) * dx; 
            particles[pt].ay += (qprefact + prefact) * dy; 
            particles[pt].az += (qprefact + prefact) * dz; 
#else
            particles[pt].ax += prefact*dx; 
            particles[pt].ay += prefact*dy; 
            particles[pt].az += prefact*dz; 
#endif
        }
    } else { // It's a leaf node
        if (node->remote == 0 && node->pt == pt) return;
        double _r = sqrt(r2 + softening2);
        double prefact = -G/(_r*_r*_r)*node->m;
        particles[pt].ax += prefact*dx; 
        particles[pt].ay += prefact*dy; 
        particles[pt].az += prefact*dz; 
    }
}

