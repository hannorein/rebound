/**
 * @file 	problem.c
 * @brief 	Example problem: Radiation forces
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example provides an implementation of the 
 * Poynting-Robertson effect. The code is using the IAS15 integrator
 * which is ideally suited for this velocity dependent force.
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "main.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "problem.h"
#include "integrator.h"
#include "integrator_ias15.h"

void additional_forces();
double betaparticles = 0.01; 	// beta parameter
				// defined as the ratio of radiation pressure over gravity

void problem_init(int argc, char* argv[]){
	// setup constants
	dt 				= 1e-3;			// initial timestep
	integrator			= IAS15;
	integrator_ias15_epsilon 	= 1e-4;			// accuracy parameter
	boxsize 			= 10;	
	tmax				= 1e5;
	N_active			= 1; 			// the star is the only massive particle
	problem_additional_forces	= additional_forces;	// setup callback function for velocity dependent forces
	init_box();
	
	
	// star is at rest at origin
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.m  = 1.;
	particles_add(star);

	// dust particles are initially on a circular orbit
	while(N<2){
		struct particle p; 
		p.m  = 0;					// massless
		double a = 1.;					// a = 1 AU
		double v = sqrt(G*(star.m*(1.-betaparticles))/a);
		double phi = tools_uniform(0,2.*M_PI);		// random phase
		p.x  = a*sin(phi);  p.y  = a*cos(phi); p.z  = 0; 
		p.vx = -v*cos(phi); p.vy = v*sin(phi); p.vz = 0;
		particles_add(p); 
	}

	system("rm -v radius.txt");					// remove previous output
}

void force_radiation(){
	const struct particle star = particles[0];				// cache
#pragma omp parallel for
	for (int i=0;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		if (p.m!=0.) continue; 						// only dust particles feel radiation forces
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 		// distance relative to star
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double c 		= 1.006491504759635e+04; 		// speed of light in unit of G=1, M_sun=1, 1year=1
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		const double F_r 	= betaparticles*G*star.m/(pr*pr);

		// Equation (5) of Burns, Lamy, Soter (1979)
		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
	}
}

void additional_forces(){
	force_radiation();							// PR drag (see above)
}

void problem_output(){
	if(output_check(400.)){						// print some information to screen
		output_timing();
	}
	if(output_check(M_PI*2.*1000.)){ 					// output radial distance every 1000 years
		FILE* f = fopen("radius.txt","a");
		const struct particle star = particles[0];
		for (int i=1;i<N;i++){
			const struct particle p = particles[i]; 
			const double prx  = p.x-star.x;
			const double pry  = p.y-star.y;
			const double prz  = p.z-star.z;
			const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
			fprintf(f,"%e\t%e\n",t,pr);
		}
		fclose(f);
	}
}

void problem_finish(){
}
