/**
 * @file 	problem.c
 * @brief 	Example problem: J2 precession
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This example presents an implementation of the J2
 * gravitational moment. The equation of motions are integrated with
 * the 15th order IAS15 integrator. The parameters in this examples 
 * have been chosen to represent those of Saturn, but you can easily
 * change them or even include higher order terms in the multipole 
 * expansion.
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
#include "input.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"
#include "problem.h"

void force_J2();
const double J2planet			= 16298e-6; 		// J2 of Saturn (Murray and Dermott p 531) 
const double Mplanet			= 0.00028588598; 	// mass of Saturn in solar masses 
const double Rplanet			= 0.00038925688; 	// radius of Saturn in AU

double ObliquityPlanet;						// obliquity of the planet

void problem_init(int argc, char* argv[]){
	// Setup constants
	integrator			= IAS15;
	dt 				= 1e-6;			// initial timestep
	boxsize 			= 0.01;	
	tmax				= 3e1;
	N_active			= 2; 			// only the star and the planet are massive.
	problem_additional_forces 	= force_J2;
	init_box();
	
	
	// Planet
	struct particle planet;
	planet.m  = Mplanet;
	planet.x  = 0; planet.y  = 0; planet.z  = 0;
	planet.vx = 0; planet.vy = 0; planet.vz = 0;
	particles_add(planet);
	// Read obliquity from command line. Default is 0. 
	ObliquityPlanet 		= input_get_double(argc,argv,"Obliquity",0.)/180.*M_PI;
	
	

	// Add one dust particles
	while(N<2){ 						// two particles in total (planet and dust particle) 
		struct particle p; 
		p.m  = 0;					// massless
		double a = Rplanet*3.;				// small distance from planet (makes J2 important)
		double e = 0.1;
		double v = sqrt((1.+e)/(1.-e)*G*planet.m/a);	// setup eccentric orbit (ignores J2)
		p.x  = (1.-e)*a; p.y  = 0; p.z  = 0; 
		p.vx = 0; p.vy = v; p.vz = 0;
		p.x += planet.x; 	p.y += planet.y; 	p.z += planet.z;
		p.vx += planet.vx; 	p.vy += planet.vy; 	p.vz += planet.vz;
		particles_add(p); 
	}
	
	tools_move_to_center_of_momentum();

	system("rm -v a.txt");					// delete previous output
}

void force_J2(){
	if (J2planet==0) return;
	// Star 
	const struct particle planet = particles[0];		// cache
#pragma omp parallel for
	for (int i=1;i<N;i++){
		const struct particle p = particles[i]; 	// cache
		const double sprx = p.x-planet.x;
		const double spry = p.y-planet.y;
		const double sprz = p.z-planet.z;
		const double prx  = sprx*cos(-ObliquityPlanet) + sprz*sin(-ObliquityPlanet);
		const double pry  = spry;
		const double prz  =-sprx*sin(-ObliquityPlanet) + sprz*cos(-ObliquityPlanet);
		const double pr2  = prx*prx + pry*pry + prz*prz; 		// distance^2 relative to planet
		const double fac  = 3.*G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);

		const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
		const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
		const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);
		
		particles[i].ax += pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet);
		particles[i].ay += pay;
		particles[i].az +=-pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet);
	}
}

void problem_output(){
	if(output_check(4000.*dt)){				// output something to screen	
		output_timing();
	}
	if(output_check(M_PI*2.*0.01)){				// output semimajor axis to file
		FILE* f = fopen("a.txt","a");
		const struct particle planet = particles[0];
		for (int i=1;i<N;i++){
			struct orbit o = tools_p2orbit(particles[i],planet);
			fprintf(f,"%.15e\t%.15e\t%.15e\t%.15e\n",t,o.a,o.e,o.omega);
		}
		fclose(f);
	}
}

void problem_finish(){
}
