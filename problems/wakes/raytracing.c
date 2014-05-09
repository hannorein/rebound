/**
 * @file 	raytracing.c
 * @brief 	Tree-based ray tracing algorithm.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein
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
#ifdef OPENMP
#include <omp.h>
#endif
#include "main.h"
#include "particle.h"
#include "output.h"
#include "boundaries.h"
#include "tree.h"
#include "tools.h"
#include "collisions.h"

int tree_does_ray_hit_particle_above(struct cell* c, double a[3], double b[3],int pr);
void tree_does_ray_hit_particle(struct cell* c, double a[3], double b[3], double sun_b[3], double* flux, double* height);

void tree_raytrace(int N_rays, double B, double phi, double sun_B, double sun_phi, double* flux, double* opacity ){
	double b[3]; // incoming ray vector
	b[0] 		= cos(phi)*cos(B);
	b[1] 		= sin(phi)*cos(B);
	b[2] 		= sin(B);
	double sun_b[3]; // light ray vector
	sun_b[0] 	= cos(sun_phi)*cos(sun_B);
	sun_b[1] 	= sin(sun_phi)*cos(sun_B);
	sun_b[2] 	= sin(sun_B);
	const double taninv = 1./tan(B);
	double raylength = taninv * boxsize_x/2.+2.*collisions_max_r; 
	int nghostray=0;
	double boxlengthray = 0;
	while(boxlengthray<raylength){
		if (nghostray==0){
			boxlengthray += boxsize_x/2.;
		}else{
			boxlengthray += boxsize_x;
		}
		nghostray++;
	}


	int N_ray_hit_particle 	= 0;	// Number of rays that hit a particle
	double F_reflected 	= 0.;	// Integrated flux. 
#pragma omp parallel for 
	for(int j=0;j<N_rays;j++){
		double a[3]; // ray position in xy plane
		a[0] = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
		a[1] = tools_uniform(-boxsize_y/2.,boxsize_y/2.); // Only really works in square boxes.
		a[2] = 0;
		int transparency = 0;
		double flux = -1;	// -1 means ray passes through ring, 0 means it hits a dark particle, 1 means it hits sees the most iluminated part of a a particle
		double height = -boxsize_z;
		for (int gbx=-nghostray; gbx<=nghostray; gbx++){
			for (int gby=-nghostray; gby<=nghostray; gby++){
				struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,0);
				double as[3]; // ishifted ray position in xy plane
				as[0] = a[0]+gb.shiftx;
				as[1] = a[1]+gb.shifty;
				as[2] = a[2];
				for (int i=0;i<root_n;i++){
					struct cell* c = tree_root[i];
					tree_does_ray_hit_particle(c,as,b, sun_b, &flux, &height);	
				}
			}
		}
		if (flux>=0.){
			N_ray_hit_particle 	+= 1;
			F_reflected 		+= flux;
		}
	}
	

	*opacity	= (double)N_ray_hit_particle	/(double)N_rays;
	*flux		= (double)F_reflected		/(double)N_rays;
}

int tree_does_ray_hit_particle_above(struct cell* c, double a[3], double b[3],int pr){
	double t1 = b[1]*(a[2]-c->z)-b[2]*(a[1]-c->y);
	double t2 = b[2]*(a[0]-c->x)-b[0]*(a[2]-c->z);
	double t3 = b[0]*(a[1]-c->y)-b[1]*(a[0]-c->x);
	double bot2 = b[0]*b[0]+b[1]*b[1]+b[2]*b[2];
	double distance2 = (t1*t1+t2*t2+t3*t3)/bot2;
	double width=c->w*0.86602540378+collisions_max_r;
	if (distance2<width*width){
		if (c->pt<0){
			for(int i=0;i<8;i++){
				struct cell* o = c->oct[i];
				if(o){
					if(tree_does_ray_hit_particle_above(o,a,b,pr)){
						return 1;
					}
				}
			}
		}else{
			struct particle p = particles[c->pt];
			double u = b[0]*(a[0]-p.x) + b[1]*(a[1]-p.y) + b[2]*(a[2]-p.z);
			double v = (p.x-a[0])*(p.x-a[0]) + (p.y-a[1])*(p.y-a[1]) + (p.z-a[2])*(p.z-a[2]);
			double s = u*u - v + p.r*p.r;
			if (s>0){ // line intersects sphere
				double d  = -u + sqrt(s); // other intersection at -u-sqrt(s). Choosing the one with larger z value.
				//double cx = a[0] + d*b[0];	// intersection point 
				//double cy = a[1] + d*b[1];
				double cz = a[2] + d*b[2];
				if (cz>a[2] && c->pt != pr){
					return 1; 	// ray hits other sphere above scattering sphere
				}
			}
		}
	}
	return 0;
}

void tree_does_ray_hit_particle(struct cell* c, double a[3], double b[3], double sun_b[3], double* flux, double* height){
	double t1 = b[1]*(a[2]-c->z)-b[2]*(a[1]-c->y);
	double t2 = b[2]*(a[0]-c->x)-b[0]*(a[2]-c->z);
	double t3 = b[0]*(a[1]-c->y)-b[1]*(a[0]-c->x);
	double bot2 = b[0]*b[0]+b[1]*b[1]+b[2]*b[2];
	double distance2 = (t1*t1+t2*t2+t3*t3)/bot2;
	double width=c->w*0.86602540378+collisions_max_r;
	if (distance2<width*width){
		if (c->pt<0){
			for(int i=0;i<8;i++){
				struct cell* o = c->oct[i];
				if(o){
					tree_does_ray_hit_particle(o,a,b, sun_b, flux, height);
				}
			}
		}else{
			struct particle p = particles[c->pt];
			double u = b[0]*(a[0]-p.x) + b[1]*(a[1]-p.y) + b[2]*(a[2]-p.z);
			double v = (p.x-a[0])*(p.x-a[0]) + (p.y-a[1])*(p.y-a[1]) + (p.z-a[2])*(p.z-a[2]);
			double s = u*u - v + p.r*p.r;
			if (s>0){ // line intersects sphere
				double d  = -u + sqrt(s); 	// There's another intersection at -u-sqrt(s). Choosing the one with larger z value.
				double cx = a[0] + d*b[0];	// intersection point 
				double cy = a[1] + d*b[1];
				double cz = a[2] + d*b[2];

				if (cz>(*height)){			// intersection happens before previous intersection
					double nx = cx - p.x; 		// normal to sphere
					double ny = cy - p.y;
					double nz = cz - p.z;
					double n  = sqrt(nx*nx + ny*ny + nz*nz);
					nx /= n;
					ny /= n;
					nz /= n;
				
					double cosNB = nx*sun_b[0] + ny*sun_b[1] + nz*sun_b[2];
					if (cosNB>0){ // dayside of particle
						// Check for shadow.
						int nghostray = 1; // check only one ghost box for now. TODO
						int isinshadow = 0;
						for (int gbx=-nghostray; gbx<=nghostray; gbx++){
							for (int gby=-nghostray; gby<=nghostray; gby++){
								struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,0);
								double inters[3]; // ray position in xy plane
								inters[0] = cx+gb.shiftx;
								inters[1] = cy+gb.shifty;
								inters[2] = cz; 
								for (int i=0;i<root_n;i++){
									struct cell* ce = tree_root[i];
									if (tree_does_ray_hit_particle_above(ce,inters,sun_b,c->pt)){
										isinshadow = 1;
									}
								}
							}
						}
						if (isinshadow){
							(*flux) 	= 0;
						}else{
							(*flux) 	= cosNB;
						}
						(*height)	= cz;
					}else{
						(*flux) 	= 0;
						(*height)	= cz;
					}
				}
			}
		}
	}
}
