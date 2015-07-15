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
#include "rebound.h"

void force_J2(struct reb_simulation* r);
const double J2planet			= 16298e-6; 		// J2 of Saturn (Murray and Dermott p 531) 
const double Mplanet			= 0.00028588598; 	// mass of Saturn in solar masses 
const double Rplanet			= 0.00038925688; 	// radius of Saturn in AU

double ObliquityPlanet			= 0.;			// obliquity of the planet

double tmax				= 3e1;			// Maximum integration time
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->integrator			= RB_IT_IAS15;
	r->dt 				= 1e-6;			// initial timestep
	r->N_active			= 2; 			// only the star and the planet are massive.
	r->additional_forces 		= force_J2;
	
	// Planet
	struct reb_particle planet;
	planet.m  = Mplanet;
	planet.x  = 0; planet.y  = 0; planet.z  = 0;
	planet.vx = 0; planet.vy = 0; planet.vz = 0;
	reb_add(r, planet);

	struct reb_particle p; 
	p.m  = 0;					// massless
	double a = Rplanet*3.;				// small distance from planet (makes J2 important)
	double e = 0.1;
	double v = sqrt((1.+e)/(1.-e)*r->G*planet.m/a);	// setup eccentric orbit (ignores J2)
	p.x  = (1.-e)*a; p.y  = 0; p.z  = 0; 
	p.vx = 0; p.vy = v; p.vz = 0;
	p.x += planet.x; 	p.y += planet.y; 	p.z += planet.z;
	p.vx += planet.vx; 	p.vy += planet.vy; 	p.vz += planet.vz;
	reb_add(r, p); 
	
	reb_move_to_com(r);

	system("rm -v a.txt");					// delete previous output

	r->heartbeat = heartbeat;
	reb_integrate(r, tmax);
}

void force_J2(struct reb_simulation* r){
	if (J2planet==0) return;
	// Star 
	const struct reb_particle planet = r->particles[0];		// cache
	const int N = r->N;
#pragma omp parallel for
	for (int i=1;i<N;i++){
		const struct reb_particle p = r->particles[i]; 	// cache
		const double sprx = p.x-planet.x;
		const double spry = p.y-planet.y;
		const double sprz = p.z-planet.z;
		const double prx  = sprx*cos(-ObliquityPlanet) + sprz*sin(-ObliquityPlanet);
		const double pry  = spry;
		const double prz  =-sprx*sin(-ObliquityPlanet) + sprz*cos(-ObliquityPlanet);
		const double pr2  = prx*prx + pry*pry + prz*prz; 		// distance^2 relative to planet
		const double fac  = 3.*r->G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);

		const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
		const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
		const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);
		
		r->particles[i].ax += pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet);
		r->particles[i].ay += pay;
		r->particles[i].az +=-pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet);
	}
}

void heartbeat(struct reb_simulation* r){
	if(reb_output_check(r, 4000.*r->dt)){				// output something to screen	
		reb_output_timing(r, tmax);
	}
	if(reb_output_check(r,M_PI*2.*0.01)){				// output semimajor axis to file
		FILE* f = fopen("a.txt","a");
		const struct reb_particle planet = r->particles[0];
		const int N = r->N;
		for (int i=1;i<N;i++){
			struct reb_orbit o = reb_tools_p2orbit(r->G, r->particles[i],planet);
			fprintf(f,"%.15e\t%.15e\t%.15e\t%.15e\n",r->t,o.a,o.e,o.omega);
		}
		fclose(f);
	}
}

