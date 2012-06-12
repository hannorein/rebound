/**
 * @file 	integrator.c
 * @brief 	RADAU15 integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the radau15 integration scheme.  
 * This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * See Everhart, 1985, ASSL Vol. 115, IAU Colloq. 83, Dynamics of 
 * Comets, Their Origin and Evolution, 185.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Dave Spiegel.
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
#ifdef TREE
#error RADAU15 integrator not working with TREE module.
#endif
#ifdef MPI
#error RADAU15 integrator not working with MPI.
#endif

int integrator_radau_init_done = 0;
double integrator_radau_accuracy = 1e-4;
void integrator_radau_init();
void integrator_radau_step();

void integrator_part1(){
	// Do nothing here. This is for the first drift part in a leapfrog-like DKD integrator.
}

// This function updates the acceleration on all particles. 
// It uses the current position and velocity data in the (struct particle*) particles structure.
// Note: this does currently not work with MPI or any TREE module.
void integrator_update_acceleration(){
  gravity_calculate_acceleration();
  if (problem_additional_forces) problem_additional_forces();
}

void integrator_part2(){
	if (!integrator_radau_init_done){
		integrator_radau_init();
		integrator_radau_init_done = 1;
	}

	integrator_radau_step();
}
double h[8], xc[8],vc[7], r[28],c[21],d[21],s[9];


unsigned int niter 	= 6;
int N3 			= 0; 	// This is just N*3
int N3allocated 	= 0; 	// Allocated memory size

double* x  = NULL;
double* v  = NULL;
double* a  = NULL;
double* x1  = NULL;
double* v1  = NULL;
double* a1  = NULL;

double* g[7];
double* b[7];
double* e[7];

struct particle* frame_out = NULL;
struct particle* frame_in  = NULL;

void integrator_radau_init() {
	// Gauss Radau spacings
	h[0] = 0.0;
	h[1] = 0.05626256053692215;
	h[2] = 0.18024069173689236;
	h[3] = 0.35262471711316964;
	h[4] = 0.54715362633055538;
	h[5] = 0.73421017721541053;
	h[6] = 0.88532094683909577;
	h[7] = 0.97752061356128750;

	//  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
	xc[0] = 0.5;
	xc[1] = 0.16666666666666667;
	xc[2] = 0.08333333333333333;
	xc[3] = 0.05;
	xc[4] = 0.03333333333333333;
	xc[5] = 0.02380952380952381;
	xc[6] = 0.01785714285714286;
	xc[7] = 0.01388888888888889;

	//  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
	vc[0] = 0.5;
	vc[1] = 0.3333333333333333;
	vc[2] = 0.25;
	vc[3] = 0.2;
	vc[4] = 0.1666666666666667;
	vc[5] = 0.1428571428571429;
	vc[6] = 0.125;

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

	for (int j=0;j<7;++j) {
		g[j] = NULL;
		b[j] = NULL;
		e[j] = NULL;
	}
}
  
void integrator_radau_realloc_memory(){
	for (int l=0;l<7;++l) {
		g[l] = realloc(g[l],sizeof(double)*N3);
		b[l] = realloc(b[l],sizeof(double)*N3);
		e[l] = realloc(e[l],sizeof(double)*N3);
		for (int k=0;k<N3;k++){
			b[l][k] = 0;
			e[l][k] = 0;
		}

	}
	//
	x = realloc(x,sizeof(double)*N3);
	v = realloc(v,sizeof(double)*N3);
	a = realloc(a,sizeof(double)*N3);
	//
	x1 = realloc(x1,sizeof(double)*N3);
	v1 = realloc(v1,sizeof(double)*N3);
	a1 = realloc(a1,sizeof(double)*N3);
	//
	frame_out = realloc(frame_out,sizeof(struct particle)*N);
	frame_in  = realloc(frame_in, sizeof(struct particle)*N);
}
  
void integrator_radau_step() {
	N3 = 3*N;
	if (N3 > N3allocated) {
		integrator_radau_realloc_memory();
		N3allocated = N3;
		niter = 6;
	}
	
	// frame_out = frame_in;
	memcpy(frame_in ,particles,N*sizeof(struct particle));
	memcpy(frame_out,particles,N*sizeof(struct particle));
	struct particle* _frame_orig = particles;



	particles = frame_out;
	integrator_update_acceleration();
	for(int k=0;k<N;++k) {
		x1[3*k]   = frame_in[k].x;
		x1[3*k+1] = frame_in[k].y;
		x1[3*k+2] = frame_in[k].z;
		//
		v1[3*k]   = frame_in[k].vx;
		v1[3*k+1] = frame_in[k].vy;
		v1[3*k+2] = frame_in[k].vz;
		//
		a1[3*k]   = frame_out[k].ax;
		a1[3*k+1] = frame_out[k].ay;  
		a1[3*k+2] = frame_out[k].az;
	}

	for(int k=0;k<N3;++k) {
		g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
		g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
		g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
		g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
		g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
		g[5][k] = b[6][k]*d[20] + b[5][k];
		g[6][k] = b[6][k];
	}

	for (int main_loop_counter=0;main_loop_counter<niter;++main_loop_counter) {
		for(int j=1;j<8;++j) {

			s[0] = dt * h[j];
			s[1] = s[0] * s[0] * 0.5;
			s[2] = s[1] * h[j] * 0.3333333333333333;
			s[3] = s[2] * h[j] * 0.5;
			s[4] = s[3] * h[j] * 0.6;
			s[5] = s[4] * h[j] * 0.6666666666666667;
			s[6] = s[5] * h[j] * 0.7142857142857143;
			s[7] = s[6] * h[j] * 0.75;
			s[8] = s[7] * h[j] * 0.7777777777777778;

			for(int k=0;k<N3;++k) {
				x[k] = 	s[8]*b[6][k] + s[7]*b[5][k] + s[6]*b[4][k] + s[5]*b[3][k] + s[4]*b[2][k] + s[3]*b[1][k] + s[2]*b[0][k] + s[1]*a1[k] + s[0]*v1[k] + x1[k];
			}
			// This should be made optional
#define VELOCITYDEPENDENDFORCE
#ifdef VELOCITYDEPENDENDFORCE
			s[0] = dt * h[j];
			s[1] = s[0] * h[j] * 0.5;
			s[2] = s[1] * h[j] * 0.6666666666666667;
			s[3] = s[2] * h[j] * 0.75;
			s[4] = s[3] * h[j] * 0.8;
			s[5] = s[4] * h[j] * 0.8333333333333333;
			s[6] = s[5] * h[j] * 0.8571428571428571;
			s[7] = s[6] * h[j] * 0.875;

			for(int k=0;k<N3;++k) {
				v[k] =  s[7]*b[6][k] + s[6]*b[5][k] + s[5]*b[4][k] + s[4]*b[3][k] + s[3]*b[2][k] + s[2]*b[1][k] + s[1]*b[0][k] + s[0]*a1[k] + v1[k];
			}
#endif // VELOCITYDEPENDENDFORCE

			for(int k=0;k<N;++k) {
				frame_out[k] = frame_in[k];
				double rr[3];
				double vv[3];
				double drr[3];
				double dvv[3];

				rr[0] = x[3*k];
				rr[1] = x[3*k+1];
				rr[2] = x[3*k+2];

				drr[0] = rr[0] - frame_in[k].x;
				drr[1] = rr[1] - frame_in[k].y;
				drr[2] = rr[2] - frame_in[k].z;
				frame_out[k].x += drr[0];
				frame_out[k].y += drr[1];
				frame_out[k].z += drr[2];

				vv[0] = v[3*k];
				vv[1] = v[3*k+1];
				vv[2] = v[3*k+2];

				dvv[0] = vv[0] - frame_in[k].vx;
				dvv[1] = vv[1] - frame_in[k].vy;
				dvv[2] = vv[2] - frame_in[k].vz;
				frame_out[k].vx += dvv[0];
				frame_out[k].vy += dvv[1];
				frame_out[k].vz += dvv[2];
			}

			/*
			HR: I do not understadn this.:

			if (interaction->IsSkippingJPLPlanets()) {
				frame_out.SetTime(frame_in+dt*h[j]);
				frame_out.ForceJPLEphemerisData();
			}
			//
			*/
			particles = frame_out;
			integrator_update_acceleration();

			for(int k=0;k<N;++k) {
				a[3*k]   = frame_out[k].ax;
				a[3*k+1] = frame_out[k].ay;  
				a[3*k+2] = frame_out[k].az;
			}
			switch (j) {
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
					} break;
				default:
					printf("Error\n");
			}

		}
	}
	double dt_done = dt;

	// Estimate suitable sequence size for the next call
	double tmp = 0.0;
	for(int k=0;k<N3;++k) {
		double _fabsb6k = fabs(b[6][k]);
		if (_fabsb6k>tmp) tmp = _fabsb6k;
	}
	// if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(dt),7));
	// if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(dt.GetDouble()),7));
	if (tmp!=0.0) tmp /= (72.0 * pow(fabs(dt),7));

	if (tmp < 1.0e-50) { // is equal to zero?
		// dt = dt_done * 1.4;
		dt = dt_done * 1.4;
	} else {

		// old rule...
		// dt = copysign(secure_pow(accuracy/tmp,0.1111111111111111),dt_done); // 1/9=0.111...
		// dt = copysign(secure_pow(accuracy/tmp,0.1111111111111111),dt_done.GetDouble()); // 1/9=0.111...
		dt = copysign(pow(integrator_radau_accuracy/tmp,0.1111111111111111),dt_done); // 1/9=0.111...

	}
	//
	// if (fabs(dt/dt_done) < 1.0) {
	if (fabs(dt/dt_done) < 1.0) {
		dt = dt_done * 0.8;
		printf("Radau: step rejected! New proposed dt: %f\n",dt);
		// HR TODO COPY!!
		particles = _frame_orig;
		memcpy(particles,frame_in,N*sizeof(struct particle));
		niter = 6;
		return;
	}
	//
	if (fabs(dt/dt_done) > 1.4) dt = dt_done * 1.4;

	// std::cerr << "RA15: new dt: " << dt.GetDouble() << std::endl;

	// Find new position and velocity values at end of the sequence
	tmp = dt_done * dt_done;
	for(int k=0;k<N3;++k) {
		x1[k] = ( xc[7]*b[6][k] +
				xc[6]*b[5][k] + 
				xc[5]*b[4][k] + 
				xc[4]*b[3][k] + 
				xc[3]*b[2][k] + 
				xc[2]*b[1][k] + 
				xc[1]*b[0][k] + 
				xc[0]*a1[k]   ) * tmp + v1[k]*dt_done + x1[k];

		v1[k] = ( vc[6]*b[6][k] + 
				vc[5]*b[5][k] + 
				vc[4]*b[4][k] +
				vc[3]*b[3][k] + 
				vc[2]*b[2][k] + 
				vc[1]*b[1][k] +
				vc[0]*b[0][k] + 
				a1[k])        * dt_done + v1[k];
	}

	for(int k=0;k<N;++k) {
		double rr[3];
		double vv[3];
		double drr[3];
		double dvv[3];

		frame_out[k] = frame_in[k];

		rr[0] = x1[3*k];
		rr[1] = x1[3*k+1];
		rr[2] = x1[3*k+2];

		drr[0] = rr[0] - frame_in[k].x;  
		drr[1] = rr[1] - frame_in[k].y;  
		drr[2] = rr[2] - frame_in[k].z;  
		frame_out[k].x += drr[0];
		frame_out[k].y += drr[1];
		frame_out[k].z += drr[2];

		vv[0] = v1[3*k];
		vv[1] = v1[3*k+1];
		vv[2] = v1[3*k+2];

		dvv[0] = vv[0] - frame_in[k].vx;
		dvv[1] = vv[1] - frame_in[k].vy;
		dvv[2] = vv[2] - frame_in[k].vz;
		frame_out[k].vx += dvv[0];
		frame_out[k].vy += dvv[1];
		frame_out[k].vz += dvv[2];
	}

	t += dt_done;
	niter = 2;
	particles = _frame_orig;
	memcpy(particles,frame_out,N*sizeof(struct particle));

	// Predict new B values to use at the start of the next sequence. The predicted
	// values from the last call are saved as E. The correction, BD, between the
	// actual and predicted values of B is applied in advance as a correction.
	//
	double q1 = dt / dt_done;
	double q2 = q1 * q1;
	double q3 = q1 * q2;
	double q4 = q2 * q2;
	double q5 = q2 * q3;
	double q6 = q3 * q3;
	double q7 = q3 * q4;

	for(int k=0;k<N3;++k) {

		s[0] = b[0][k] - e[0][k];
		s[1] = b[1][k] - e[1][k];
		s[2] = b[2][k] - e[2][k];
		s[3] = b[3][k] - e[3][k];
		s[4] = b[4][k] - e[4][k];
		s[5] = b[5][k] - e[5][k];
		s[6] = b[6][k] - e[6][k];

		// Estimate B values for the next sequence

		e[0][k] = q1*(b[6][k]* 7.0 + b[5][k]* 6.0 + b[4][k]* 5.0 + b[3][k]* 4.0 + b[2][k]* 3.0 + b[1][k]*2.0 + b[0][k]);
		e[1][k] = q2*(b[6][k]*21.0 + b[5][k]*15.0 + b[4][k]*10.0 + b[3][k]* 6.0 + b[2][k]* 3.0 + b[1][k]);
		e[2][k] = q3*(b[6][k]*35.0 + b[5][k]*20.0 + b[4][k]*10.0 + b[3][k]* 4.0 + b[2][k]);
		e[3][k] = q4*(b[6][k]*35.0 + b[5][k]*15.0 + b[4][k]* 5.0 + b[3][k]);
		e[4][k] = q5*(b[6][k]*21.0 + b[5][k]* 6.0 + b[4][k]);
		e[5][k] = q6*(b[6][k]* 7.0 + b[5][k]);
		e[6][k] = q7* b[6][k];

		b[0][k] = e[0][k] + s[0];
		b[1][k] = e[1][k] + s[1];
		b[2][k] = e[2][k] + s[2];
		b[3][k] = e[3][k] + s[3];
		b[4][k] = e[4][k] + s[4];
		b[5][k] = e[5][k] + s[5];
		b[6][k] = e[6][k] + s[6];

	}

}
