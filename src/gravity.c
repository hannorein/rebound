/**
 * @file 	gravity.c
 * @brief 	Direct gravity calculation, O(N^2).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This is the crudest implementation of an N-body code
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
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "rebound.h"
#include "tree.h"
#include "boundary.h"

#ifdef MPI
#include "communication_mpi.h"
#endif

/**
  * @brief The function loops over all trees to call calculate_forces_for_particle_from_cell() tree to calculate forces for each particle.
  * @param r REBOUND simulation to consider
  * @param pt Index of the particle the force is calculated for.
  * @param gb Ghostbox plus position of the particle (precalculated). 
  */
static void reb_calculate_acceleration_for_particle(const struct reb_simulation* const r, const int pt, const struct reb_ghostbox gb);

/**
 * Main Gravity Routine
 */
void reb_calculate_acceleration(struct reb_simulation* r){
	struct reb_particle* const particles = r->particles;
	const int N = r->N;
	const int N_active = r->N_active;
	const double G = r->G;
	const double softening2 = r->softening*r->softening;
	const unsigned int _gravity_ignore_terms = r->gravity_ignore_terms;
	const int _N_active = ((N_active==-1)?N:N_active) - r->N_var;
	const int _N_real   = N  - r->N_var;
	const int _testparticle_type   = r->testparticle_type;
	switch (r->gravity){
		case REB_GRAVITY_NONE: // Do nothing.
		break;
		case REB_GRAVITY_BASIC:
		{
			const int nghostx = r->nghostx;
			const int nghosty = r->nghosty;
			const int nghostz = r->nghostz;
#pragma omp parallel for schedule(guided)
			for (int i=0; i<N; i++){
				particles[i].ax = 0; 
				particles[i].ay = 0; 
				particles[i].az = 0; 
			}
			// Summing over all Ghost Boxes
			for (int gbx=-nghostx; gbx<=nghostx; gbx++){
			for (int gby=-nghosty; gby<=nghosty; gby++){
			for (int gbz=-nghostz; gbz<=nghostz; gbz++){
				struct reb_ghostbox gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
				// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
				for (int i=0; i<_N_real; i++){
				for (int j=0; j<_N_active; j++){
					if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0) )) continue;
					if (_gravity_ignore_terms==2 && ((j==0 || i==0) )) continue;
					if (i==j) continue;
					const double dx = (gb.shiftx+particles[i].x) - particles[j].x;
					const double dy = (gb.shifty+particles[i].y) - particles[j].y;
					const double dz = (gb.shiftz+particles[i].z) - particles[j].z;
					const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
					const double prefact = -G/(_r*_r*_r)*particles[j].m;
					
					particles[i].ax    += prefact*dx;
					particles[i].ay    += prefact*dy;
					particles[i].az    += prefact*dz;
				}
				}
                if (_testparticle_type){
				for (int i=0; i<_N_active; i++){
				for (int j=_N_active; j<_N_real; j++){
					if (_gravity_ignore_terms==1 && ((j==1 && i==0) )) continue;
					if (_gravity_ignore_terms==2 && ((j==0 || i==0) )) continue;
					const double dx = (gb.shiftx+particles[i].x) - particles[j].x;
					const double dy = (gb.shifty+particles[i].y) - particles[j].y;
					const double dz = (gb.shiftz+particles[i].z) - particles[j].z;
					const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
					const double prefact = -G/(_r*_r*_r)*particles[j].m;
					
					particles[i].ax    += prefact*dx;
					particles[i].ay    += prefact*dy;
					particles[i].az    += prefact*dz;
				}
				}
                }
			}
			}
			}
		}
		break;
		case REB_GRAVITY_COMPENSATED:
		{
			if (r->gravity_cs_allocatedN<N){
				r->gravity_cs = realloc(r->gravity_cs,N*sizeof(struct reb_vec3d));
				r->gravity_cs_allocatedN = N;
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
			for (int gbx=-r->nghostx; gbx<=r->nghostx; gbx++){
			for (int gby=-r->nghosty; gby<=r->nghosty; gby++){
			for (int gbz=-r->nghostz; gbz<=r->nghostz; gbz++){
				// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
				for (int i=0; i<N; i++){
					struct reb_ghostbox gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
					// Precalculated shifted position
					gb.shiftx += particles[i].x;
					gb.shifty += particles[i].y;
					gb.shiftz += particles[i].z;
					reb_calculate_acceleration_for_particle(r, i, gb);
				}
			}
			}
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
	const int N = r->N;
	const int _N_real   = N - r->N_var;
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
		case REB_GRAVITY_BASIC:
            for (int v=0;v<r->var_config_N;v++){
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
                        for (int i=0; i<_N_real; i++){
                        for (int j=i+1; j<_N_real; j++){
                            if (_gravity_ignore_terms==1 && ((j==1 && i==0) || (i==1 && j==0))) continue;
                            if (_gravity_ignore_terms==2 && ((j==0 || i==0))) continue;
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
                            dax 	+=        3.* r5inv * dk2dx * rdk1
                                    + 3.* r5inv * dk1dx * rdk2
                                    + 3.* r5inv    * dx * dk1dk2  
                                        - 15.      * dx * r7inv * rdk1 * rdk2;
                            day 	+=        3.* r5inv * dk2dy * rdk1
                                    + 3.* r5inv * dk1dy * rdk2
                                    + 3.* r5inv    * dy * dk1dk2  
                                        - 15.      * dy * r7inv * rdk1 * rdk2;
                            daz 	+=        3.* r5inv * dk2dz * rdk1
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
                            dax 	+=        3.* r5inv * dk2dx * rdk1
                                    + 3.* r5inv * dk1dx * rdk2
                                    + 3.* r5inv    * dx * dk1dk2  
                                        - 15.      * dx * r7inv * rdk1 * rdk2;
                            day 	+=        3.* r5inv * dk2dy * rdk1
                                    + 3.* r5inv * dk1dy * rdk2
                                    + 3.* r5inv    * dy * dk1dk2  
                                        - 15.      * dy * r7inv * rdk1 * rdk2;
                            daz 	+=        3.* r5inv * dk2dz * rdk1
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
static void reb_calculate_acceleration_for_particle_from_cell(const struct reb_simulation* const r, const int pt, const struct reb_treecell *node, const struct reb_ghostbox gb);

static void reb_calculate_acceleration_for_particle(const struct reb_simulation* const r, const int pt, const struct reb_ghostbox gb) {
	for(int i=0;i<r->root_n;i++){
		struct reb_treecell* node = r->tree_root[i];
		if (node!=NULL){
			reb_calculate_acceleration_for_particle_from_cell(r, pt, node, gb);
		}
	}
}

static void reb_calculate_acceleration_for_particle_from_cell(const struct reb_simulation* r, const int pt, const struct reb_treecell *node, const struct reb_ghostbox gb) {
	const double G = r->G;
	const double softening2 = r->softening*r->softening;
	struct reb_particle* const particles = r->particles;
	const double dx = gb.shiftx - node->mx;
	const double dy = gb.shifty - node->my;
	const double dz = gb.shiftz - node->mz;
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
			double mrr 	= dx*dx*node->mxx 	+ dy*dy*node->myy 	+ dz*dz*node->mzz
					+ 2.*dx*dy*node->mxy 	+ 2.*dx*dz*node->mxz 	+ 2.*dy*dz*node->myz; 
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
		if (node->pt == pt) return;
		double _r = sqrt(r2 + softening2);
		double prefact = -G/(_r*_r*_r)*node->m;
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
	}
}

