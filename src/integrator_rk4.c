/**
 * @file 	integrator.c
 * @brief 	Euler integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the Euler integration scheme.  
 * There is practically no use for this first order integrator. For 
 * non-rotating coordinates systems integrator_leapfrog should be used. 
 * For shearing sheet boundary conditions integrator_sei is well suited.
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
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "problem.h"
#include "output.h"
#include "tools.h"



// These variables have no effect for euler.
int integrator_force_is_velocitydependent 	= 1;
double integrator_epsilon 			= 0;
extern double integrator_min_dt 		= 0;


/**
 * This part has no function, as the euler scheme does not have any sub-
 * timestep.
 */
void integrator_part1(){
}

// This function updates the acceleration on all particles. 
// It uses the current position and velocity data in the (struct particle*) particles structure.
// Note: this does currently not work with MPI or any TREE module.
void integrator_update_acceleration(){
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)
	PROFILING_START()
	gravity_calculate_acceleration();
	if (problem_additional_forces) problem_additional_forces();
	PROFILING_STOP(PROFILING_CAT_GRAVITY)
	PROFILING_START()
}


double rk4(double(*f)(double, double), double dx, double x, double y)
{
	double	k1 = dx * f(x, y),
		k2 = dx * f(x + dx / 2, y + k1 / 2),
		k3 = dx * f(x + dx / 2, y + k2 / 2),
		k4 = dx * f(x + dx, y + k3);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}
 
double rate(double x, double y)
{
	return x * sqrt(y);
}




void integrator_part2(){
    double* x0  	= NULL;			// Temporary buffer for position (used for initial values at h=0) 
    double* v0  	= NULL;			//                      velocity
    double* v1  	= NULL;			//                      velocity
    double* v2  	= NULL;			//                      velocity
    double* v3  	= NULL;			//                      velocity
    double* a0  	= NULL;			//                      acceleration
    double* a1  	= NULL;			//                      acceleration
    double* a2  	= NULL;			//                      acceleration
    double* a3  	= NULL;			//                      acceleration
    
	const int N3 = 3*N;
    x0 = realloc(x0,sizeof(double)*N3);
    v0 = realloc(v0,sizeof(double)*N3);
    v1 = realloc(v1,sizeof(double)*N3);
    v2 = realloc(v2,sizeof(double)*N3);
    v3 = realloc(v3,sizeof(double)*N3);
    a0 = realloc(a0,sizeof(double)*N3);
    a1 = realloc(a1,sizeof(double)*N3);
    a2 = realloc(a2,sizeof(double)*N3);
    a3 = realloc(a3,sizeof(double)*N3);
    
    // Save initial state
	for(int i=0;i<N;i++) {
		x0[3*i]   = particles[i].x;
		x0[3*i+1] = particles[i].y;
		x0[3*i+2] = particles[i].z;
		v0[3*i]   = particles[i].vx;
		v0[3*i+1] = particles[i].vy;
		v0[3*i+2] = particles[i].vz;
		a0[3*i]   = particles[i].ax;
		a0[3*i+1] = particles[i].ay;  
		a0[3*i+2] = particles[i].az;
	}
    
    // K2
	for (int i=0;i<N;i++){
		particles[i].x  = x0[3*i+0] + 0.5* dt * v0[3*i+0];
		particles[i].y  = x0[3*i+1] + 0.5* dt * v0[3*i+1];
		particles[i].z  = x0[3*i+2] + 0.5* dt * v0[3*i+2];
		particles[i].vx = v0[3*i+0] + 0.5* dt * a0[3*i+0];
		particles[i].vy = v0[3*i+1] + 0.5* dt * a0[3*i+1];
		particles[i].vz = v0[3*i+2] + 0.5* dt * a0[3*i+2];
	}
	t += 0.5 * dt;
    integrator_update_acceleration();				// Calculate forces at interval n
	t -= 0.5 * dt;
    // Save velocities and accelerations
	for(int i=0;i<N;i++) {
		v1[3*i]   = particles[i].vx;
		v1[3*i+1] = particles[i].vy;
		v1[3*i+2] = particles[i].vz;
		a1[3*i]   = particles[i].ax;
		a1[3*i+1] = particles[i].ay;  
		a1[3*i+2] = particles[i].az;
	}
    
    

    // Reset particles position
    for(int i=0;i<N;++i) {
        particles[i].x = x0[3*i+0];	// Set inital position
        particles[i].y = x0[3*i+1];
        particles[i].z = x0[3*i+2];
    }
    
    // K3
	for (int i=0;i<N;i++){
		particles[i].x  = x0[3*i+0] + 0.5* dt * v1[3*i+0];
		particles[i].y  = x0[3*i+1] + 0.5* dt * v1[3*i+1];
		particles[i].z  = x0[3*i+2] + 0.5* dt * v1[3*i+2];
		particles[i].vx = v0[3*i+0] + 0.5* dt * a1[3*i+0];
		particles[i].vy = v0[3*i+1] + 0.5* dt * a1[3*i+1];
		particles[i].vz = v0[3*i+2] + 0.5* dt * a1[3*i+2];
	}
	t += 0.5 * dt;
    integrator_update_acceleration();				// Calculate forces at interval n
	t -= 0.5 * dt;
    // Save velocities and accelerations
	for(int i=0;i<N;i++) {
		v2[3*i]   = particles[i].vx;
		v2[3*i+1] = particles[i].vy;
		v2[3*i+2] = particles[i].vz;
		a2[3*i]   = particles[i].ax;
		a2[3*i+1] = particles[i].ay;  
		a2[3*i+2] = particles[i].az;
	}
    


    // Reset particles position
    for(int i=0;i<N;++i) {
        particles[i].x = x0[3*i+0];	// Set inital position
        particles[i].y = x0[3*i+1];
        particles[i].z = x0[3*i+2];
    }
    

    // K4 
	for (int i=0;i<N;i++){
		particles[i].x  = x0[3*i+0] + dt * v2[3*i+0];
		particles[i].y  = x0[3*i+1] + dt * v2[3*i+1];
		particles[i].z  = x0[3*i+2] + dt * v2[3*i+2];
		particles[i].vx = v0[3*i+0] + dt * a2[3*i+0];
		particles[i].vy = v0[3*i+1] + dt * a2[3*i+1];
		particles[i].vz = v0[3*i+2] + dt * a2[3*i+2];
	}
        
    t += dt;
    integrator_update_acceleration();				// Calculate forces at interval n
    // Save accelerations
	for(int i=0;i<N;i++) {
		v3[3*i]   = particles[i].vx;
		v3[3*i+1] = particles[i].vy;
		v3[3*i+2] = particles[i].vz;
		a3[3*i]   = particles[i].ax;
		a3[3*i+1] = particles[i].ay;  
		a3[3*i+2] = particles[i].az;
	}
    
    
    // Calculate
    for(int i=0;i<N;++i) {
        particles[i].x = x0[3*i+0] + dt * ( 1./6. * (v0[3*i+0] + 2.*(v1[3*i+0] + v2[3*i+0]) + v3[3*i+0]) );
        particles[i].y = x0[3*i+1] + dt * ( 1./6. * (v0[3*i+1] + 2.*(v1[3*i+1] + v2[3*i+1]) + v3[3*i+1]) );
        particles[i].z = x0[3*i+2] + dt * ( 1./6. * (v0[3*i+2] + 2.*(v1[3*i+2] + v2[3*i+2]) + v3[3*i+2]) );
        particles[i].vx = v0[3*i+0] + dt * ( 1./6. * (a0[3*i+0] + 2.*(a1[3*i+0] + a2[3*i+0]) + a3[3*i+0]) );
        particles[i].vy = v0[3*i+1] + dt * ( 1./6. * (a0[3*i+1] + 2.*(a1[3*i+1] + a2[3*i+1]) + a3[3*i+1]) );
        particles[i].vz = v0[3*i+2] + dt * ( 1./6. * (a0[3*i+2] + 2.*(a1[3*i+2] + a2[3*i+2]) + a3[3*i+2]) );
    }

}
	

