/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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

// These variables have no effect for leapfrog.
int integrator_force_is_velocitydependent 	= 1;
double integrator_epsilon 			= 0;
double integrator_min_dt 			= 0;

// Fast factorial lookup table
static const double factorial[] = {1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880., 3628800., 39916800., 479001600., 6227020800., 87178291200., 1307674368000., 20922789888000., 355687428096000., 6402373705728000., 121645100408832000., 2432902008176640000., 51090942171709440000., 1124000727777607680000., 25852016738884976640000., 620448401733239439360000., 15511210043330985984000000., 403291461126605635584000000., 10888869450418352160768000000., 304888344611713860501504000000., 8841761993739701954543616000000., 265252859812191058636308480000000., 8222838654177922817725562880000000., 263130836933693530167218012160000000., 8683317618811886495518194401280000000., 295232799039604140847618609643520000000.};

double c_n_series(unsigned int n, double z){
	double c_n = 0.;
	for (unsigned int j=0;j<13;j++){
		double term = pow(-z,j)/factorial[n+2*j];
		c_n += term;
		if (fabs(term/c_n) < 1e-17) break; // Stop if new term smaller than machine precision
	}
	return c_n;
}

double c(unsigned int n, double z){
	if (z>0.5){
		switch(n){
			case 0:
			{
				double cn4 = c(3,z/4.)*(1.+c(1,z/4.))/8.;
				double cn2 = 1./2.-z*cn4;
				double cn0 = 1.-z*cn2;
				return cn0;
			}
			case 1:
			{
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4))/16.;
				double cn3 = 1./6.-z*cn5;
				double cn1 = 1.-z*cn3;
				return cn1;
			}
			case 2:
			{
				double cn4 = c(3,z/4.)*(1.+c(1,z/4.))/8.;
				double cn2 = 1./2.-z*cn4;
				return cn2;
			}
			case 3:
			{
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4))/16.;
				double cn3 = 1./6.-z*cn5;
				return cn3;
			}
			case 4:
			{
				double cn4 = c(3,z/4.)*(1.+c(1,z/4.))/8.;
				return cn4;
			}
			case 5:
			{
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4))/16.;
				return cn5;
			}
		}
	}
	return c_n_series(n,z);
}

double integrator_G(unsigned int n, double beta, double X){
	return pow(X,n)*c(n,beta*X*X);
}

void kepler_step(){
	struct particle p0 = particles[0];
	struct particle p1 = particles[1];

	double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
	double v2 =  p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
	double beta = 2.*p0.m/r0 - v2;
	double eta = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;
	double zeta = p0.m - beta*r0;
	

	// Find solution for X using bisection method. Not good because:
	// 1) slow
	// 2) only work if period = 2pi (as of now)
	// 3) probably doesn't work for hyperbolic orbits
	// 4) might not work for retrograde orbits
	
	double X_min = 0;
	double X_max = 2.*M_PI;
	while (X_max-X_min>1e-17){
		double X = (X_max + X_min)/2.;
		double s= r0*X + eta*integrator_G(2,beta,X) + zeta*integrator_G(3,beta,X)-dt;
		if (s<0.){
			X_min = X;
		}else{
			X_max = X;
		}
	}
	double X = X_min;

	FILE* f = fopen("tmp.txt","w");
	fclose(f);
	exit(0);

	
}


// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void integrator_part1(){
	kepler_step();
	printf("\n C(t)=%f",c_n_series(2,t));
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	t+=dt/2.;
}
void integrator_part2(){
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		particles[i].vx += dt * particles[i].ax;
		particles[i].vy += dt * particles[i].ay;
		particles[i].vz += dt * particles[i].az;
		particles[i].x  += 0.5* dt * particles[i].vx;
		particles[i].y  += 0.5* dt * particles[i].vy;
		particles[i].z  += 0.5* dt * particles[i].vz;
	}
	t+=dt/2.;
}
	

