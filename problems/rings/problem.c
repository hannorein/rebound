/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate particles on a circular orbit in a fixed 
 * potential.
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "main.h"
#include "problem.h"
#include "input.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

// Star and planet (note, those wont be updated after they have been inserted)
struct particle star;
struct particle planet;
extern double 	integrator_accuracy;	// Desired accuracy. Play with this, make sure you get a converged results.
extern double 	integrator_min_dt;	// Timestep floor.
extern int 	integrator_adaptive_timestep;
void additional_forces();
const double   Rjup   = 0.00046732617; // Radius of Jupiter in AU
const double   J2jup  = 14736e-6; 	// J2 of Jupiter 

const double Rplanet = Rjup;
const double J2planet = J2jup;
int Nmoon = 0;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 				= 1e-2;
	//integrator_accuracy 		= 1e-4;
	integrator_adaptive_timestep	= 0;
	boxsize 			= 12.;
	N_active      			= 2;
	problem_additional_forces 	= additional_forces;
	init_box();

	// Initial conditions
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1;
	particles_add(star);

	// Planet 1 at 1 AU
	planet.m  = 1.e-3;
	planet.x  = 5.2; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	output_double("planet mass", planet.m);
	output_double("planet a", planet.x);

	// Moon 1 Amalthea
	struct particle moon1 = planet;
	moon1.m  = 0;
	double r 	= 0.0012122499*.5;
	double v	= sqrt(G*planet.m / r);
	double inclination = 0.37/180.*M_PI;
	moon1.x  +=  r;
	moon1.vy +=  v*cos(inclination);
	moon1.vz +=  v*sin(inclination);
	particles_add(moon1);
	Nmoon++;
	
	// Moon 2 Thebe
	struct particle moon2 = planet;
	moon2.m  = 0;
	r 	  =  0.0014833099*.5;
	v	  = sqrt(G*planet.m / r);
	inclination = 1.09/180.*M_PI;
	moon2.x  +=  r;
	moon2.vy +=  v*cos(inclination);
	moon2.vz +=  v*sin(inclination);
	particles_add(moon2);
	Nmoon++;

	system("cat config.log");
	tools_move_to_center_of_momentum();
}
		
void force_J2(){
	// Star 
	const struct particle planet = particles[1];				// cache
#pragma omp parallel for
	for (int i=2;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double prx  = p.x-planet.x;
		const double pry  = p.y-planet.y;
		const double prz  = p.z-planet.z;
		const double pr2   = prx*prx + pry*pry + prz*prz; 	// distance^2 relative to planet
		const double fac  = 3.*G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);
		particles[i].ax += fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
		particles[i].ay += fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
		particles[i].az += fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);
	}
}


void force_radiation(){
	// Star 
	const struct particle star = particles[0];				// cache
//#define SHADOW
#ifdef SHADOW
	const struct particle planet = particles[1];				// cache
#endif // SHADOW
#pragma omp parallel for
	for (int i=2+Nmoon;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star

#ifdef SHADOW
		const double plrx = planet.x-star.x;
		const double plry = planet.y-star.y;
		const double plrz = planet.z-star.z;
		const double plr  = sqrt(plrx*plrx + plry*plry + plrz*plrz); 	// distance relative to star
		int IN_SHADOW = 0;
		// Find out if particle is in the shadow of the planet
		// i.e., find out if the radial vector to the particle passes through the planet
		// If the particle isn't as far as the planet's center, it's definitely not in shadow.
		if (pr < plr)
			IN_SHADOW = 0;
		else {
			// If the angle between the particle's position vector relative to the star and the center of
			// the planet's position vector relative to the center of the star is bigger than Rplanet/pr, then
			// the planet is not in shadow.  Otherwise, it is.
			const double max_angle_from_planet_center = Rjup / pr; // planet radius / distance from star
			// Call vector to planet center "Pl" and vector to particle "p"
			// angle from planet center is arccos[Pl dot p / (plr * pr)]
			const double dotprod = prx*plrx + pry*plry + prz*plrz;
			const double cos_angle_from_planet_center = dotprod / (plr*pr);
			const double angle_from_planet_center = acos(cos_angle_from_planet_center);
			if (angle_from_planet_center < max_angle_from_planet_center)
				IN_SHADOW = 1;
		}
		if (IN_SHADOW == 1) continue;
#endif // SHADOW

		// No shadow; calculate force modification:
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double beta	= 0.012;    				 // appropriate for 100-micron dust
		const double c 		= 10064.915; 				// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		const double F_r 	= beta*G*star.m/(pr*pr);

		// Old: Equation (5) of Burns, Lamy, Soter (1979)
		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);

		// New: TODO: Why this is not the same as above. 
		// Set up tangential velocity
		// vtan + vr = v so vtan = v - vr
		/*
		const double vtan_x = prvx - (rdot * (prx/pr));
		const double vtan_y = prvy - (rdot * (pry/pr));
		const double vtan_z = prvz - (rdot * (prz/pr));

		const double modify_mass_term_x = F_r * (prx/pr); // This is the radiation pressure term that scales
		const double modify_mass_term_y = F_r * (pry/pr); // with (1/r^2) and therefore just makes the effective
		const double modify_mass_term_z = F_r * (prz/pr); // mass of the star slightly less.

		const double doppler_term_x = modify_mass_term_x * (-rdot/c);
		const double doppler_term_y = modify_mass_term_y * (-rdot/c);
		const double doppler_term_z = modify_mass_term_z * (-rdot/c);

		const double tangential_term_x = -F_r * (vtan_x/c);
		const double tangential_term_y = -F_r * (vtan_y/c);
		const double tangential_term_z = -F_r * (vtan_z/c);

		const double radiation_force_x = modify_mass_term_x + doppler_term_x + tangential_term_x;
		const double radiation_force_y = modify_mass_term_y + doppler_term_y + tangential_term_y;
		const double radiation_force_z = modify_mass_term_z + doppler_term_z + tangential_term_z;
		particles[i].ax += radiation_force_x;
		particles[i].ay += radiation_force_y;
		particles[i].az += radiation_force_z;
		*/
	}
}

void additional_forces(){
	force_J2();
//	force_radiation();
}
						
void output_a(){
	FILE* of = fopen("a.txt","w"); 
	struct particle planet = particles[1];
	for (int i=2+Nmoon;i<N;i++){
		struct particle p = particles[i];
		const double prx  = p.x-planet.x;
		const double pry  = p.y-planet.y;
		const double prz  = p.z-planet.z;
		const double pr2   = prx*prx + pry*pry + prz*prz; 	// distance^2 relative to planet
		fprintf(of,"%e\n",sqrt(pr2));
	}
	fclose(of);
}

void problem_output(){
	// Ring particles
	const int Nringparticlemax = 1000;
	if (N-2-Nmoon < Nringparticlemax) {
		printf("adding particls\n");
		struct particle ringparticle;
		double vesc;

		// Amalthea
		ringparticle = particles[2];
		vesc = 0.002954419;
		ringparticle.vx += vesc * tools_normal(1.);
		ringparticle.vy += vesc * tools_normal(1.);
		ringparticle.vz += vesc * tools_normal(1.);
		particles_add(ringparticle);

		// Thebe
		ringparticle = particles[3];
		vesc = 0.001745793;
		ringparticle.vx += vesc * tools_normal(1.);
		ringparticle.vy += vesc * tools_normal(1.);
		ringparticle.vz += vesc * tools_normal(1.);
		particles_add(ringparticle);
	}
	
	
	// Output stuff
	if(output_check(2.*M_PI)){
		output_timing();
	}
	if(output_check(.5)){
		output_a();
	}
}
void problem_inloop(){
}

void problem_finish(){
}
