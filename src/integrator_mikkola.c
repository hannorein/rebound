/**
 * @file 	integrator.c
 * @brief 	Mikkola integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the Mikkola integration scheme.  
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein
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

// Fast inverse factorial lookup table
static const double invfactorial[] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

double c_n_series(unsigned int n, double z){
	double c_n = 0.;
	for (unsigned int j=0;j<13;j++){
		double term = pow(-z,j)*invfactorial[n+2*j];
		c_n += term;
		if (fabs(term/c_n) < 1e-17) break; // Stop if new term smaller than machine precision
	}
	return c_n;
}

double c(unsigned int n, double z){
	if (z>0.5){
		// Speed up conversion with 4-folding formula
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
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4.))/16.;
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
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4.))/16.;
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
				double cn5 = (c(5,z/4.)+c(4,z/4.)+c(3,z/4.)*c(2,z/4.))/16.;
				return cn5;
			}
		}
	}
	return c_n_series(n,z);
}

double integrator_G(unsigned int n, double beta, double X){
	return pow(X,n)*c(n,beta*X*X);
}

void kepler_step(int i){
	double M = particles[0].m;
	struct particle p1 = particles[i];

	double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
	double v2 =  p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
	double beta = 2.*M/r0 - v2;
	double eta = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;
	double zeta = M - beta*r0;
	

	double X = 0;  // TODO: find a better initial estimate.
	double G1,G2,G3;
	for (int n_hg=0;n_hg<20;n_hg++){
		G2 = integrator_G(2,beta,X);
		G3 = integrator_G(3,beta,X);
		G1 = X-beta*G3;
		double s   = r0*X + eta*G2 + zeta*G3-dt;
		double sp  = r0 + eta*G1 + zeta*G2;
		double dX  = -s/sp; // Newton's method
		
		//double G0 = 1.-beta*G2;
		//double spp = r0 + eta*G0 + zeta*G1;
		//double dX  = -(s*sp)/(sp*sp-0.5*s*spp); // Householder 2nd order formula
		X+=dX;
		if (fabs(dX/X)<1e-15) break; // TODO: Make sure loop does not get stuck (add maximum number of iterations) 
	}

	double r = r0 + eta*G1 + zeta*G2;
	double f = 1.-M*G2/r0;
	double g = dt - M*G3;
	double fd = -M*G1/(r0*r); 
	double gd = 1.-M*G2/r; 

	particles[i].x = f*p1.x + g*p1.vx;
	particles[i].y = f*p1.y + g*p1.vy;
	particles[i].z = f*p1.z + g*p1.vz;

	particles[i].vx = fd*p1.x + gd*p1.vx;
	particles[i].vy = fd*p1.y + gd*p1.vy;
	particles[i].vz = fd*p1.z + gd*p1.vz;

}


void integrator_part1(){
}
void integrator_to_jacobi(){
	double M = 0.;
	double capx  = 0.; 
	double capy  = 0.; 
	double capz  = 0.; 
	double capvx = 0.; 
	double capvy = 0.; 
	double capvz = 0.; 

	for (int i=0;i<N;i++){
		double m = particles[i].m;
		M += m;
		if (i>0){
			particles[i].x  -= capx;
			particles[i].y  -= capy;
			particles[i].z  -= capz;
			particles[i].vx -= capvx;
			particles[i].vy -= capvy;
			particles[i].vz -= capvz;
		}
		capx  += particles[i].x* m/M;
		capy  += particles[i].y* m/M;
		capz  += particles[i].z* m/M;
		capvx += particles[i].vx*m/M;
		capvy += particles[i].vy*m/M;
		capvz += particles[i].vz*m/M;
	}
}
void integrator_to_heliocentric(){
	double M = 0.;
	double capx  = 0.; 
	double capy  = 0.; 
	double capz  = 0.; 
	double capvx = 0.; 
	double capvy = 0.; 
	double capvz = 0.; 
	for (int i=0;i<N;i++){
		double m = particles[i].m;
		M += m;
		capx  += particles[i].x* m/M;  
		capy  += particles[i].y* m/M;
		capz  += particles[i].z* m/M;
		capvx += particles[i].vx*m/M;
		capvy += particles[i].vy*m/M;
		capvz += particles[i].vz*m/M;
	}
	for (int i=N-1;i>0;i--){
		double m = particles[i].m;
		capx  -= particles[i].x* m/M; 
		capy  -= particles[i].y* m/M;
		capz  -= particles[i].z* m/M;
		capvx -= particles[i].vx*m/M;
		capvy -= particles[i].vy*m/M;
		capvz -= particles[i].vz*m/M;
		M -= m;
			
		particles[i].x  += capx;
		particles[i].y  += capy;
		particles[i].z  += capz;
		particles[i].vx += capvx;
		particles[i].vy += capvy;
		particles[i].vz += capvz;
	}
}

void integrator_part2(){
	integrator_to_jacobi();
	integrator_to_heliocentric();
	for (int i=1;i<N;i++){
		kepler_step(i);
	}
	t+=dt;
}
	

