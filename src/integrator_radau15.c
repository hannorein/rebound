/**
 * @file 	integrator.c
 * @brief 	RADAU15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the radau15 integration scheme.  
 * This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * See Everhart, 1985, ASSL Vol. 115, IAU Colloq. 83, Dynamics of 
 * Comets, Their Origin and Evolution, 185.
 * Part of this code is based a function from the ORSE package.
 * See orsa.sourceforge.net for more details on their implementation.
 *
 * The user might want to change the following variables in the 
 * problem.c file:
 * extern int integrator_force_is_velocitydependent;
 * extern double integrator_epsilon;
 * extern double integrator_min_dt;
 * 
 * 
 * @section 	LICENSE
 * Copyright (c) 2011-2012 Hanno Rein, Dave Spiegel.
 * Copyright (c) 2002-2004 Pasquale Tricarico.
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
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "problem.h"
#include "output.h"
#include "tools.h"

#ifdef TREE
#error RADAU15 integrator not working with TREE module.
#endif
#ifdef MPI
#error RADAU15 integrator not working with MPI.
#endif

int 	integrator_force_is_velocitydependent	= 1;	// Turn this off to safe some time if the force is not velocity dependent.
double 	integrator_epsilon 			= 0;	// Magnitude of last term in series expansion devided by the acceleration is smaller than this value or timestep is rejected. 
							// Play with integrator_epsilon to make sure you get a converged results. 
							// The true fractional error is often many orders of magnitude smaller.
							// If it is zero, then a constant timestep is used (default). 
double 	integrator_min_dt 			= 0;	// Minimum timestep used as a floor when adaptive timestepping is enabled.
double	integrator_error			= 0;	// Error estimate in last timestep (used for debugging only)
unsigned int integrator_iterations_max		= 10;	// Maximum number of iterations in predictor/corrector loop
unsigned long integrator_iterations_max_exceeded= 0;	// Count how many times the iteration did not converge


const double h[8]	= { 0.0, 0.05626256053692215, 0.18024069173689236, 0.35262471711316964, 0.54715362633055538, 0.73421017721541053, 0.88532094683909577, 0.97752061356128750}; // Gauss Radau spacings
const double xc[8] 	= { 0.5, 0.16666666666666667, 0.08333333333333333, 0.05, 0.03333333333333333, 0.02380952380952381, 0.01785714285714286, 0.01388888888888889}; // 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
const double vc[7] 	= { 0.5, 0.3333333333333333, 0.25, 0.2, 0.1666666666666667, 0.1428571428571429, 0.125}; // 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8

double r[28],c[21],d[21],s[9]; // These constants will be set dynamically.

int N3allocated 		= 0; 	// Size of allocated arrays.
int integrator_radau_init_done 	= 0;	// Calculate coefficients once.

double* x   = NULL;	// Temporary buffer for position
double* v   = NULL;	//                      velocity
double* a   = NULL;	//                      acceleration
double* x1  = NULL;	// Temporary buffer for position
double* v1  = NULL;	//                      velocity
double* a1  = NULL;	//                      acceleration

double* g[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* b[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* e[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;

struct particle* particles_out = NULL; // Temporary particle buffer.

void integrator_part1(){
	// Do nothing here. This is only used in a leapfrog-like DKD integrator.
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

int integrator_radau_step(); // Does the actual timestep.

void integrator_part2(){
	if (!integrator_radau_init_done){ 	// Generate coefficients.
		int l=0;
		for (int j=1;j<8;++j) {
			for(int k=0;k<j;++k) {
				r[l] = 1.0 / (h[j] - h[k]);
				++l;
			}
		}
		c[0] = -h[1];
		d[0] =  h[1];
		l=0;
		for (int j=2;j<7;++j) {
			++l;
			c[l] = -h[j] * c[l-j+1];
			d[l] =  h[1] * d[l-j+1];
			for(int k=2;k<j;++k) {
				++l;
				c[l] = c[l-j] - h[j] * c[l-j+1];
				d[l] = d[l-j] + h[k] * d[l-j+1];
			}
			++l;
			c[l] = c[l-j] - h[j];
			d[l] = d[l-j] + h[j]; 
		}
		integrator_radau_init_done = 1;
	}
	// Try until a step was successful.
	while(!integrator_radau_step());
}
 
int integrator_radau_step() {
	const int N3 = 3*N;
	if (N3 > N3allocated) {
		for (int l=0;l<7;++l) {
			g[l] = realloc(g[l],sizeof(double)*N3);
			b[l] = realloc(b[l],sizeof(double)*N3);
			e[l] = realloc(e[l],sizeof(double)*N3);
			for (int k=0;k<N3;k++){
				b[l][k] = 0;
				e[l][k] = 0;
			}
		}
		x = realloc(x,sizeof(double)*N3);
		v = realloc(v,sizeof(double)*N3);
		a = realloc(a,sizeof(double)*N3);
		x1 = realloc(x1,sizeof(double)*N3);
		v1 = realloc(v1,sizeof(double)*N3);
		a1 = realloc(a1,sizeof(double)*N3);
		particles_out = realloc(particles_out,sizeof(struct particle)*N);
		N3allocated = N3;
	}
	
	struct particle* particles_in  = particles;
	// integrator_update_acceleration(); // Not needed. Forces are already calculated in main routine.

	for(int k=0;k<N;k++) {
		x1[3*k]   = particles[k].x;
		x1[3*k+1] = particles[k].y;
		x1[3*k+2] = particles[k].z;
		v1[3*k]   = particles[k].vx;
		v1[3*k+1] = particles[k].vy;
		v1[3*k+2] = particles[k].vz;
		a1[3*k]   = particles[k].ax;
		a1[3*k+1] = particles[k].ay;  
		a1[3*k+2] = particles[k].az;
	}

	for(int k=0;k<N3;k++) {
		g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
		g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
		g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
		g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
		g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
		g[5][k] = b[6][k]*d[20] + b[5][k];
		g[6][k] = b[6][k];
	}

	double predictor_corrector_error = 1;
	int iterations = 0;
	while(predictor_corrector_error>1e-10){						// Predictor corrector loop
		if (iterations>=integrator_iterations_max){
			integrator_iterations_max_exceeded++;
			if (integrator_iterations_max_exceeded==1){
				fprintf(stderr,"\n\033[1mWarning!\033[0m A large number (>100) of predictor corrector loops in integrator_radau15.c did not converge. This is typically an indication of the timestep being too large.\n");
			}
			break;								// Quit predictor corrector loop
		}
		predictor_corrector_error = 0;
		iterations++;

		for(int n=1;n<8;n++) {							// Loop over interval using Gauss-Radau spacings

			s[0] = dt * h[n];
			s[1] = s[0] * s[0] * 0.5;
			s[2] = s[1] * h[n] * 0.3333333333333333;
			s[3] = s[2] * h[n] * 0.5;
			s[4] = s[3] * h[n] * 0.6;
			s[5] = s[4] * h[n] * 0.6666666666666667;
			s[6] = s[5] * h[n] * 0.7142857142857143;
			s[7] = s[6] * h[n] * 0.75;
			s[8] = s[7] * h[n] * 0.7777777777777778;

			for(int k=0;k<N3;k++) {						// Predict positions at interval n using b values
				x[k] = 	s[8]*b[6][k] + s[7]*b[5][k] + s[6]*b[4][k] + s[5]*b[3][k] + s[4]*b[2][k] + s[3]*b[1][k] + s[2]*b[0][k] + s[1]*a1[k] + s[0]*v1[k] + x1[k];
			}
			
			if (integrator_force_is_velocitydependent){
				s[0] = dt * h[n];
				s[1] = s[0] * h[n] * 0.5;
				s[2] = s[1] * h[n] * 0.6666666666666667;
				s[3] = s[2] * h[n] * 0.75;
				s[4] = s[3] * h[n] * 0.8;
				s[5] = s[4] * h[n] * 0.8333333333333333;
				s[6] = s[5] * h[n] * 0.8571428571428571;
				s[7] = s[6] * h[n] * 0.875;

				for(int k=0;k<N3;k++) {					// Predict velocities at interval n using b values
					v[k] =  s[7]*b[6][k] + s[6]*b[5][k] + s[5]*b[4][k] + s[4]*b[3][k] + s[3]*b[2][k] + s[2]*b[1][k] + s[1]*b[0][k] + s[0]*a1[k] + v1[k];
				}
			}

			// Prepare particles arrays for force calculation
			for(int k=0;k<N;k++) {
				particles_out[k] = particles_in[k];

				particles_out[k].x = x[3*k+0];
				particles_out[k].y = x[3*k+1];
				particles_out[k].z = x[3*k+2];

				particles_out[k].vx = v[3*k+0];
				particles_out[k].vy = v[3*k+1];
				particles_out[k].vz = v[3*k+2];
			}

			particles = particles_out;
			integrator_update_acceleration();				// Calculate forces at interval n

			for(int k=0;k<N;++k) {
				a[3*k]   = particles[k].ax;
				a[3*k+1] = particles[k].ay;  
				a[3*k+2] = particles[k].az;
			}
			switch (n) {							// Improve b and g values
				case 1: 
					for(int k=0;k<N3;++k) {
						double tmp = g[0][k];
						g[0][k]  = (a[k] - a1[k]) * r[0];
						b[0][k] += g[0][k] - tmp;
					} break;
				case 2: 
					for(int k=0;k<N3;++k) {
						double tmp = g[1][k];
						double gk = a[k] - a1[k];
						g[1][k] = (gk*r[1] - g[0][k])*r[2];
						tmp = g[1][k] - tmp;
						b[0][k] += tmp * c[0];
						b[1][k] += tmp;
					} break;
				case 3: 
					for(int k=0;k<N3;++k) {
						double tmp = g[2][k];
						double gk = a[k] - a1[k];
						g[2][k] = ((gk*r[3] - g[0][k])*r[4] - g[1][k])*r[5];
						tmp = g[2][k] - tmp;
						b[0][k] += tmp * c[1];
						b[1][k] += tmp * c[2];
						b[2][k] += tmp;
					} break;
				case 4:
					for(int k=0;k<N3;++k) {
						double tmp = g[3][k];
						double gk = a[k] - a1[k];
						g[3][k] = (((gk*r[6] - g[0][k])*r[7] - g[1][k])*r[8] - g[2][k])*r[9];
						tmp = g[3][k] - tmp;
						b[0][k] += tmp * c[3];
						b[1][k] += tmp * c[4];
						b[2][k] += tmp * c[5];
						b[3][k] += tmp;
					} break;
				case 5:
					for(int k=0;k<N3;++k) {
						double tmp = g[4][k];
						double gk = a[k] - a1[k];
						g[4][k] = ((((gk*r[10] - g[0][k])*r[11] - g[1][k])*r[12] - g[2][k])*r[13] - g[3][k])*r[14];
						tmp = g[4][k] - tmp;
						b[0][k] += tmp * c[6];
						b[1][k] += tmp * c[7];
						b[2][k] += tmp * c[8];
						b[3][k] += tmp * c[9];
						b[4][k] += tmp;
					} break;
				case 6:
					for(int k=0;k<N3;++k) {
						double tmp = g[5][k];
						double gk = a[k] - a1[k];
						g[5][k] = (((((gk*r[15] - g[0][k])*r[16] - g[1][k])*r[17] - g[2][k])*r[18] - g[3][k])*r[19] - g[4][k])*r[20];
						tmp = g[5][k] - tmp;
						b[0][k] += tmp * c[10];
						b[1][k] += tmp * c[11];
						b[2][k] += tmp * c[12];
						b[3][k] += tmp * c[13];
						b[4][k] += tmp * c[14];
						b[5][k] += tmp;
					} break;
				case 7:
					for(int k=0;k<N3;++k) {
						double tmp = g[6][k];
						double gk = a[k] - a1[k];
						g[6][k] = ((((((gk*r[21] - g[0][k])*r[22] - g[1][k])*r[23] - g[2][k])*r[24] - g[3][k])*r[25] - g[4][k])*r[26] - g[5][k])*r[27];
						tmp = g[6][k] - tmp;	
						b[0][k] += tmp * c[15];
						b[1][k] += tmp * c[16];
						b[2][k] += tmp * c[17];
						b[3][k] += tmp * c[18];
						b[4][k] += tmp * c[19];
						b[5][k] += tmp * c[20];
						b[6][k] += tmp;
						
						// Monitor change in b[6][k]/a[k]. The iteration is converged if it is close to 0.
						double fabstmp = fabs(tmp/a[k]);
						if (fabstmp>predictor_corrector_error && isfinite(fabstmp)){
							predictor_corrector_error = fabstmp;
						}
					} break;
			}
		//	if (n==7 && iterations>7)printf("%d %e\n",iterations,predictor_corrector_error);
		}
	}
		//	printf("\n");
	const double dt_done = dt;
	if (integrator_epsilon>0){
		// Estimate error (given by last term in series expansion) 
		double error = 0.0;
		for(int k=0;k<N3;++k) {
			double errork = fabs(b[6][k]/a[k]);
			if (!isnan(errork) && !isinf(errork) && errork>error) error = errork;
		}
		integrator_error = error; // Only used for debugging and monitoring
		// Do not change timestep if all accelerations equal to zero.
		if  (error>0.0){
			double dt_new = pow(integrator_epsilon/error,1./7.)*dt_done; 
			const double safety_factor = 0.75;  // Empirically chosen so that timestep are occasionally rejected but not too often.
			// New timestep is smaller.
			if (fabs(dt_new/dt_done) < 1.0) {
				dt_new = dt_done * safety_factor;
				if (dt_new<integrator_min_dt) dt_new = integrator_min_dt;
				if (dt_done>integrator_min_dt){
					particles = particles_in;
					dt = dt_new;
					return 0; // Step rejected. Do again. 
				}
			}else{
				// New timestep is larger.
				dt_new *= safety_factor;	// This safety_factor doesn't necesarily have to be the same as above
				if (dt_new<integrator_min_dt) dt_new = integrator_min_dt;
				if (fabs(dt_new/dt_done) > 1./safety_factor) dt_new = dt_done /safety_factor;	// Don't increase the timestep by too much compared to the last one.
			}
			dt = dt_new;
		}
	}

	// Find new position and velocity values at end of the sequence
	const double dt_done2 = dt_done * dt_done;
	for(int k=0;k<N3;++k) {
		x1[k] += (xc[7]*b[6][k] + xc[6]*b[5][k] + xc[5]*b[4][k] + xc[4]*b[3][k] + xc[3]*b[2][k] + xc[2]*b[1][k] + xc[1]*b[0][k] + xc[0]*a1[k]) 
			* dt_done2 + v1[k] * dt_done;

		v1[k] += (vc[6]*b[6][k] + vc[5]*b[5][k] + vc[4]*b[4][k] + vc[3]*b[3][k] + vc[2]*b[2][k] + vc[1]*b[1][k] + vc[0]*b[0][k] + a1[k])
			* dt_done;
	}

	t += dt_done;
	// Swap particle buffers
	particles = particles_in;

	for(int k=0;k<N;++k) {
		particles[k].x = x1[3*k+0];	// Set final position
		particles[k].y = x1[3*k+1];
		particles[k].z = x1[3*k+2];

		particles[k].vx = v1[3*k+0];	// Set final velocity
		particles[k].vy = v1[3*k+1];
		particles[k].vz = v1[3*k+2];
	}


	// Predict new B values to use at the start of the next sequence. The predicted
	// values from the last call are saved as E. The correction, BD, between the
	// actual and predicted values of B is applied in advance as a correction.
	//
	const double q1 = dt / dt_done;
	const double q2 = q1 * q1;
	const double q3 = q1 * q2;
	const double q4 = q2 * q2;
	const double q5 = q2 * q3;
	const double q6 = q3 * q3;
	const double q7 = q3 * q4;

	for(int k=0;k<N3;++k) {
		double be0 = b[0][k] - e[0][k];
		double be1 = b[1][k] - e[1][k];
		double be2 = b[2][k] - e[2][k];
		double be3 = b[3][k] - e[3][k];
		double be4 = b[4][k] - e[4][k];
		double be5 = b[5][k] - e[5][k];
		double be6 = b[6][k] - e[6][k];

		e[0][k] = q1*(b[6][k]* 7.0 + b[5][k]* 6.0 + b[4][k]* 5.0 + b[3][k]* 4.0 + b[2][k]* 3.0 + b[1][k]*2.0 + b[0][k]);
		e[1][k] = q2*(b[6][k]*21.0 + b[5][k]*15.0 + b[4][k]*10.0 + b[3][k]* 6.0 + b[2][k]* 3.0 + b[1][k]);
		e[2][k] = q3*(b[6][k]*35.0 + b[5][k]*20.0 + b[4][k]*10.0 + b[3][k]* 4.0 + b[2][k]);
		e[3][k] = q4*(b[6][k]*35.0 + b[5][k]*15.0 + b[4][k]* 5.0 + b[3][k]);
		e[4][k] = q5*(b[6][k]*21.0 + b[5][k]* 6.0 + b[4][k]);
		e[5][k] = q6*(b[6][k]* 7.0 + b[5][k]);
		e[6][k] = q7* b[6][k];

		b[0][k] = e[0][k] + be0;
		b[1][k] = e[1][k] + be1;
		b[2][k] = e[2][k] + be2;
		b[3][k] = e[3][k] + be3;
		b[4][k] = e[4][k] + be4;
		b[5][k] = e[5][k] + be5;
		b[6][k] = e[6][k] + be6;
	}
	return 1; // Success.
}
