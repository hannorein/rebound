/**
 * @file 	problem.c
 * @brief 	Example problem: Radiation forces
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example shows how to integrate circumplanetary
 * dust particles using the `integrator_ias15.c` module.
 * The example sets the function pointer `problem_additional_forces`
 * to its own function that describes the radiation forces.
 * The example uses a beta parameter of 0.01. 
 * The output is custom too, outputting the semi-major axis of 
 * every dust particle relative to the planet. 
 * Only one dust particle is used in this example, but there could be
 * many.
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
#include "integrator.h"
#include "problem.h"

void force_radiation();
double betaparticles = 0.01; 	// Beta parameter. 
				// Defined as the ratio of radiation pressure over gravity.

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 			= 1e-4;	// Initial timestep.
	boxsize 		= 10;	
	tmax			= 1e6;
	N_active		= 2; 	// Only the star and the planet are massive.
	problem_additional_forces 	= force_radiation;
	init_box();
	
	
	// Star
	struct particle star;
	star.x  = 0; star.y  = 0; star.z  = 0;
	star.vx = 0; star.vy = 0; star.vz = 0;
	star.ax = 0; star.ay = 0; star.az = 0;
	star.m  = 1.;
	particles_add(star);


	// Planet 
	struct particle planet;
	planet.m  = 1e-3;
	planet.x  = 1; planet.y  = 0.; planet.z  = 0.;
	planet.ax = 0; planet.ay = 0; planet.az = 0;
	planet.vx = 0;
	planet.vy = sqrt(G*(star.m+planet.m)/planet.x);
	planet.vz = 0;
	particles_add(planet);
	
	

	// Dust particles
	while(N<3){ 	// Three particles in total (star, planet, dust particle) 
		struct particle p; 
		p.m  = 0;		// massless
		double r = 0.001;	// distance from planet planet
		double v = sqrt(G*planet.m/r);
		p.x  = r; p.y  = 0; p.z  = 0; 
		p.vx = 0; p.vy = v; p.vz = 0;
		p.x += planet.x; 	p.y += planet.y; 	p.z += planet.z;
		p.vx += planet.vx; 	p.vy += planet.vy; 	p.vz += planet.vz;
		p.ax = 0; p.ay = 0; p.az = 0;
		particles_add(p); 
	}
	
	tools_move_to_center_of_momentum();

	system("rm -v a.txt");	
}

void force_radiation(){
	const struct particle star = particles[0];				// cache
#pragma omp parallel for
	for (int i=0;i<N;i++){
		const struct particle p = particles[i]; 			// cache
		if (p.m!=0.) continue; 						// Only dust particles feel radiation forces
		const double prx  = p.x-star.x;
		const double pry  = p.y-star.y;
		const double prz  = p.z-star.z;
		const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
		const double prvx = p.vx-star.vx;
		const double prvy = p.vy-star.vy;
		const double prvz = p.vz-star.vz;

		const double c 		= 1.006491504759635e+04; 		// speed of light.
		const double rdot 	= (prvx*prx + prvy*pry + prvz*prz)/pr; 	// radial velocity relative to star
		const double F_r 	= betaparticles*G*star.m/(pr*pr);

		// Equation (5) of Burns, Lamy, Soter (1979)
		particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
		particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
		particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
	}
}

void problem_output(){
	if(output_check(M_PI*2.)){
		output_timing();
	}
	if(output_check(M_PI*2.*100.)){ // output every 100 years
		FILE* f = fopen("a.txt","a");
		const struct particle planet = particles[1];
		for (int i=2;i<N;i++){
			const struct particle p = particles[i]; 
			const double prx  = p.x-planet.x;
			const double pry  = p.y-planet.y;
			const double prz  = p.z-planet.z;
			const double pr   = sqrt(prx*prx + pry*pry + prz*prz); 	// distance relative to star
			
			const double pvx  = p.vx-planet.vx;
			const double pvy  = p.vy-planet.vy;
			const double pvz  = p.vz-planet.vz;
			const double pv   = sqrt(pvx*pvx + pvy*pvy + pvz*pvz); 	// distance relative to star
			
			const double a = -G*planet.m/( pv*pv - 2.*G*planet.m/pr );			// semi major axis
			
			fprintf(f,"%e\t%e\n",t,a);
		}
		fclose(f);
	}
}

void problem_finish(){
}
