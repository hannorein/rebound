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
	const unsigned int _gravity_ignore_10 = r->gravity_ignore_10;
	const int _N_start  = (r->integrator==REB_INTEGRATOR_WH?1:0);
	const int _N_active = ((N_active==-1)?N:N_active) - r->N_var;
	const int _N_real   = N  - r->N_var;
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
				for (int i=_N_start; i<_N_real; i++){
				for (int j=_N_start; j<_N_active; j++){
					if (_gravity_ignore_10 && j==1 && i==0 ) continue;
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
#pragma omp parallel for schedule(guided)
			for (int i=_N_start; i<_N_active; i++){
			for (int j=i+1; j<_N_active; j++){
				if (_gravity_ignore_10 && j==1 && i==0 ) continue;
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
#pragma omp parallel for schedule(guided)
			for (int i=_N_active; i<_N_real; i++){
			for (int j=_N_start; j<_N_active; j++){
				if (_gravity_ignore_10 && i==1 && j==0 ) continue;
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
	const unsigned int _gravity_ignore_10 = r->gravity_ignore_10;
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
		// No break!
		case REB_GRAVITY_BASIC:
#pragma omp parallel for schedule(guided)
			for (int i=_N_real; i<N; i++){
				particles[i].ax = 0.; 
				particles[i].ay = 0.; 
				particles[i].az = 0.; 
			}
#pragma omp parallel for schedule(guided)
			for (int i=_N_real; i<N; i++){
			for (int j=i+1; j<N; j++){
				if (_gravity_ignore_10 && ((i==_N_real+1 && j==_N_real) || (j==_N_real+1 && i==_N_real)) ) continue;
				const double dx = particles[i-N/2].x - particles[j-N/2].x;
				const double dy = particles[i-N/2].y - particles[j-N/2].y;
				const double dz = particles[i-N/2].z - particles[j-N/2].z;
				const double r2 = dx*dx + dy*dy + dz*dz;
				const double _r  = sqrt(r2);
				const double r3inv = 1./(r2*_r);
				const double r5inv = 3.*r3inv/r2;
				const double ddx = particles[i].x - particles[j].x;
				const double ddy = particles[i].y - particles[j].y;
				const double ddz = particles[i].z - particles[j].z;
				const double Gmi = G * particles[i].m;
				const double Gmj = G * particles[j].m;
				
				const double dax =   ddx * ( dx*dx*r5inv - r3inv )
						   + ddy * ( dx*dy*r5inv )
						   + ddz * ( dx*dz*r5inv );
				const double day =   ddx * ( dy*dx*r5inv )
						   + ddy * ( dy*dy*r5inv - r3inv )
						   + ddz * ( dy*dz*r5inv );
				const double daz =   ddx * ( dz*dx*r5inv )
						   + ddy * ( dz*dy*r5inv )
						   + ddz * ( dz*dz*r5inv - r3inv );
				
				particles[i].ax += Gmj * dax;
				particles[i].ay += Gmj * day;
				particles[i].az += Gmj * daz;
				
				particles[j].ax -= Gmi * dax;
				particles[j].ay -= Gmi * day;
				particles[j].az -= Gmi * daz;
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

