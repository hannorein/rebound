/**
 * @file 	librebound.c
 * @brief 	Declarations and functions for shared library access to the IAS15 and MIKKOLA integrators.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Dave Spiegel <dave@ias.edu>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel
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
#include <sys/time.h>
#include <string.h>
#include "particle.h"
#include "gravity.h"
#include "particle.h"
#include "main.h"
#include "librebound.h"
#include "boundaries.h"
#include "integrator.h"
#include "integrator_mikkola.h"

// Default values of parameters and constants
double dt 	= 0.01;	
double t 	= 0;
double tmax	= 0;
double G 	= 1;
double softening = 0;
double timing 	= 0;
extern int Nmax;	
extern int iter;  // TODO DEBUG
const char *build_str = "Built on: " __DATE__ " " __TIME__;

// Function pointer to additional forces
void (*problem_additional_forces) () = NULL;

// Particle getter/setter methods.
void setp(struct particle* _p){
	free(particles);
	particles = malloc(sizeof(struct particle)*N);
	memcpy(particles,_p,sizeof(struct particle)*N);
}
struct particle particle_get(int i){
	return particles[i];
}
struct particle* particles_get(){
	return particles;
}
void set_additional_forces(void (* _cb)()){
	problem_additional_forces = _cb;
}

// Integrate for 1 step
void rebound_step(){ 
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	if (N<=0){
		fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
		return;
	}
	integrator_part1();
	gravity_calculate_acceleration();
	if (N_megno){
		gravity_calculate_variational_acceleration();
	}
	if (problem_additional_forces) problem_additional_forces();
	integrator_part2();
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	timing = timing_final-timing_initial;
}

void integrator_set(int i){
	integrator = i;
}

// Integrate for 1 step
void reset(){ 
	dt 		= 0.01;
	t 		= 0;
	tmax		= 0;
	G 		= 1;
	softening 	= 0;
	N 		= 0;
	Nmax 		= 0;
	N_active 	= -1;
	N_megno 	= 0;
	iter		= 0;
	timing		= 0.;
	free(particles);
	particles 	= NULL;
	integrator_reset();
	struct timeval tim;
	gettimeofday(&tim, NULL);
	srand ( tim.tv_usec + getpid());
}


// Integrate until t=_tmax (or slightly more if exactFinishTime=0)
void integrate(double _tmax, int exactFinishTime, int keepSynchronized){
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initial = tim.tv_sec+(tim.tv_usec/1000000.0);
	tmax = _tmax;
	double dt_last_done = dt;
	int last_step = 0;
	int integrator_mikkola_synchronize_manually_init = integrator_mikkola_synchronize_manually;
	int integrator_mikkola_persistent_particles_init = integrator_mikkola_persistent_particles;
	integrator_mikkola_particles_modified = 1;
	if (N_megno || keepSynchronized){
		integrator_mikkola_synchronize_manually = 0;
		integrator_mikkola_persistent_particles = 0;
	}else{
		integrator_mikkola_synchronize_manually = 1;
		integrator_mikkola_persistent_particles = 1;
	}
	if (problem_additional_forces==NULL){
		integrator_force_is_velocitydependent = 0;
	}
	while(t<tmax && last_step<2){
		if (N<=0){
			fprintf(stderr,"\n\033[1mError!\033[0m No particles found. Exiting.\n");
			return;
		}
		integrator_part1();
		gravity_calculate_acceleration();
		if (N_megno){
			gravity_calculate_variational_acceleration();
		}
		if (problem_additional_forces) problem_additional_forces();
		integrator_part2();
		
		if (t+dt>=tmax && exactFinishTime==1){
			integrator_synchronize();
			dt = tmax-t;
			last_step++;
		}else{
			dt_last_done = dt;
		}
	}
	integrator_synchronize();
	dt = dt_last_done;
	integrator_mikkola_synchronize_manually = integrator_mikkola_synchronize_manually_init;
	integrator_mikkola_persistent_particles = integrator_mikkola_persistent_particles_init;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	timing = timing_final-timing_initial;
}
	 
