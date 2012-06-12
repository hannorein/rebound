/**
 * @file 	problem.c
 * @brief 	Example problem: radau15.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the 15th order RADAU15 integrator
 * to integrate the outer planets of the solar system. The initial 
 * conditions are taken from Applegate et al 1986. Pluto is a test
 * particle.
 * 
 * @section 	LICENSE
 * Copyright (c) 2012 Hanno Rein, Shangfei Liu, Dave Spiegel
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
#include "problem.h"
#include "main.h"
#include "input.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"
#include "integrator.h"

double ss_pos[6][3] = 
{
	{-4.06428567034226e-3,	-6.08813756435987e-3,	-1.66162304225834e-6	}, // Sun
	{+3.40546614227466e+0,	+3.62978190075864e+0,	+3.42386261766577e-2	}, // Jupiter
	{+6.60801554403466e+0,	+6.38084674585064e+0,	-1.36145963724542e-1	}, // Saturn
	{+1.11636331405597e+1,	+1.60373479057256e+1,	+3.61783279369958e-1	}, // Uranus
	{-3.01777243405203e+1,	+1.91155314998064e+0,	-1.53887595621042e-1	}, // Neptune
	{-2.13858977531573e+1,	+3.20719104739886e+1,	+2.49245689556096e+0	}  // Pluto
};
double ss_vel[6][3] = 
{
	{+6.69048890636161e-6,	-6.33922479583593e-6,	-3.13202145590767e-9	}, // Sun
	{-5.59797969310664e-3,	+5.51815399480116e-3,	-2.66711392865591e-6	}, // Jupiter
	{-4.17354020307064e-3,	+3.99723751748116e-3,	+1.67206320571441e-5	}, // Saturn
	{-3.25884806151064e-3,	+2.06438412905916e-3,	-2.17699042180559e-5	}, // Uranus
	{-2.17471785045538e-4,	-3.11361111025884e-3,	+3.58344705491441e-5	}, // Neptune
	{-1.76936577252484e-3,	-2.06720938381724e-3,	+6.58091931493844e-4	}  // Pluto
};

double ss_mass[6] =
{
	1.00000597682, 	// Sun + inner planets
	1./1047.355,	// Jupiter
	1./3501.6,	// Saturn
	1./22869.,	// Uranus
	1./19314.,	// Neptune
	0.		// Pluto
};

const double k	 	= 0.01720209895;	// Gaussian constant 
#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

#ifdef INTEGRATOR_RADAU15
extern int integrator_adaptive_timestep;	
extern int integrator_force_is_velocitydependend;
extern double integrator_accuracy;
#endif // INTEGRATOR_RADAU15

void velocity_dependend_force();

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",40);	// Inital timestep in days
	N_active	= 5;
	tmax		= 7.3e10;		// 200 Myr
	G		= k*k;
	problem_additional_forces = velocity_dependend_force;
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(200); 			// Init box with width 200 astronomical units

	// Initial conditions
	for (int i=0;i<6;i++){
		struct particle p;
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = ss_mass[i];
		particles_add(p); 
	}
#ifdef INTEGRATOR_RADAU15
	integrator_adaptive_timestep 		= 1;	
	integrator_force_is_velocitydependend	= 0;
	integrator_accuracy			= 1e-5;
#endif // INTEGRATOR_RADAU15
#ifdef INTEGRATOR_WH
	// Move to heliocentric frame (required by WHM)
	for (int i=1;i<N;i++){
		particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
		particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
	}
	particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
	particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
#endif // INTEGRATOR_WH

	system("rm -f orbits.txt");
}

void velocity_dependend_force(){
	//RADAU15 can handle arbitrary velocity dependend forces such as these:
	//for (int i=1;i<N;i++){
	//	particles[i].ax -= 0.0000001*particles[i].vx;
	//	particles[i].ay -= 0.0000001*particles[i].vy;
	//	particles[i].az -= 0.0000001*particles[i].vz;
	//}
}

void problem_inloop(){
}

// Calculate the relative energy error as a diagnostic tool.
double energy_init = 0;
void output_energy(){
	double energy_kin = 0;
	double energy_grav = 0; 
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		energy_kin += 1./2.*p.m*(p.vx*p.vx + p.vy*p.vy + p.vz*p.vz); 
		for (int j=i+1;j<N;j++){
			struct particle p2 = particles[j];
			double dx = p.x - p2.x;
			double dy = p.y - p2.y;
			double dz = p.z - p2.z;
			energy_grav += -G*p.m*p2.m/sqrt(dx*dx + dy*dy + dz*dz);
		}
	}
	if (energy_init==0){
		energy_init = energy_kin + energy_grav;
	}
	
	FILE* of = fopen("energy.txt","a+"); 
	if (of==NULL){
		printf("\n\nError while opening file.\n");
		return;
	}
	fprintf(of,"%e\t%e\t%e\n",t, fabs((energy_kin + energy_grav - energy_init)/energy_init), energy_kin+energy_grav);
	fclose(of);
}

void problem_output(){
	if (output_check(10000.)){
		output_energy();
		output_timing();
	}
	if (output_check(3652422.)){ 	// output heliocentric orbital elements every 10000 years
		output_append_orbits("orbits.txt");
	}
}

void problem_finish(){
}
