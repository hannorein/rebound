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
extern double   Rjup   = 7.1492e9; // Radius of Jupiter in cm
extern double   Rearth = 6.3781e8; // Radius of Earth in cm
void radiation_forces();

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 				= 1e-1;
	// integrator_min_dt		= 1e-2;		// Set a timestep floor.
	integrator_accuracy 		= 1e-4;
	boxsize 			= 2.8;
	N_active      			= 2;
	problem_additional_forces 	= radiation_forces;
	init_box();

	// Initial conditions
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1;
	particles_add(star);

	// Planet 1 at 1 AU
	planet.m  = 1.e-3;
	planet.r  = Rjup; // set the planet's radius
	planet.x  = 1.; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	output_double("planet mass", planet.m);
	output_double("planet a", planet.x);

	// Ring particles
	double planet_hill = planet.x*pow(planet.m/3./star.m,1./3.);
	output_double("planet hill", planet_hill);
	double r_inner =  planet_hill*input_get_double(argc, argv, "r_inner", 0.1);
	double r_outer =  planet_hill*input_get_double(argc, argv, "r_outer", 0.4);
	output_double("ring inner", r_inner);
	output_double("ring outer", r_outer);
	int _N = 1000;
	for (int i = 0; i < _N; i++) {
		struct particle ringparticle = planet;
		ringparticle.m  = 0;
		double r 	= tools_uniform(r_inner,r_outer);
		double v	= sqrt(G*planet.m / r);
		double phi	= tools_uniform(0,2.*M_PI);
		ringparticle.x  +=  r*cos(phi);
		ringparticle.y  +=  r*sin(phi);
		ringparticle.vx += -v*sin(phi);
		ringparticle.vy +=  v*cos(phi);
		particles_add(ringparticle);
	}
	system("cat config.log");
	tools_move_to_center_of_momentum();
}

void radiation_forces(){
	// Star 
	const struct particle star = particles[0];				// cache
	const struct particle planet = particles[1];				// cache
#pragma omp parallel for
	for (int i=2;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
		const double plrx = planet.x-star.x;
		const double plry = planet.y-star.y;
		const double plrz = planet.z-star.z;
		const double plr  = sqrt(plrx*plrx + plry*plry + plrz*plrz); 	// distance relative to star

		const int USE_SHADOW = 1;
		int IN_SHADOW = 0;
		if (USE_SHADOW == 1)
		  {
		    // Find out if particle is in the shadow of the planet
		    // i.e., find out if the radial vector to the particle passes through the planet
		    // If the particle isn't as far as the planet's center, it's definitely not in shadow.
		    if (pr < plr)
		      IN_SHADOW = 0;
		    else {
		      // If the angle between the particle's position vector relative to the star and the center of
		      // the planet's position vector relative to the center of the star is bigger than Rplanet/pr, then
		      // the planet is not in shadow.  Otherwise, it is.
		      const double max_angle_from_planet_center = planet.r / pr; // planet radius / distance from star
		      // Call vector to planet center "Pl" and vector to particle "p"
		      // angle from planet center is arccos[Pl dot p / (plr * pr)]
		      const dotprod = prx*plrx + pry*plry + prz*plrz;
		      const cos_angle_from_planet_center = dotprod / (plr*pr);
		      const double angle_from_planet_center = acos(cos_angle_from_planet_center);
		      if (angle_from_planet_center < max_angle_from_planet_center)
			IN_SHADOW = 1;
		    }
		  }
		if (IN_SHADOW == 1) return;

		// No shadow; calculate force modification:
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double beta	= 0.012;     // appropriate for 100-micron dust
		const double c 		= 10064.915; 				// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		// why does rdot divide by pr?
		const double F_r 	= beta*G*star.m/(pr*pr);

		// Set up tangential velocity
		// vtan + vr = v so vtan = v - vr
		const double vtan_x = prvx - (rdot * (prx/pr));
		const double vtan_y = prvy - (rdot * (pry/pr));
		const double vtan_z = prvz - (rdot * (prz/pr));

		const double modify_mass_term_x = F_r * (prx/pr); // This is the radiation pressure term that scales
		const double modify_mass_term_y = F_r * (pry/pr); // with (1/r^2) and therefore just makes the effective
		const double modify_mass_term_z = F_r * (prz/pr); // mass of the star slightly less.

		const double doppler_term_x = modify_mass_term_x * (-rdot/c);
		const double doppler_term_y = modify_mass_term_y * (-rdot/c);
		const double doppler_term_z = modify_mass_term_z * (-rdot/c);

		const double tangential_term_x = F_r * (vtan_x/c);
		const double tangential_term_y = F_r * (vtan_y/c);
		const double tangential_term_z = F_r * (vtan_z/c);

		const double radiation_force_x = modify_mass_term_x + doppler_term_x + tangential_term_x;
		const double radiation_force_y = modify_mass_term_y + doppler_term_y + tangential_term_y;
		const double radiation_force_z = modify_mass_term_z + doppler_term_z + tangential_term_z;
		// Old:
		//particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		//particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		//particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);

		// New:
		particles[i].ax += radiation_force_x;
		particles[i].ay += radiation_force_y;
		particles[i].az += radiation_force_z;
	}
}


void problem_output(){
	if(output_check(2.*M_PI)){
		output_timing();
	}
	if(output_check(124.)){
		//output_append_ascii("position.txt");
		//output_append_orbits("orbits.txt");
	}
}
void problem_inloop(){
}

void problem_finish(){
}
