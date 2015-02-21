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

struct particle* p_j  = NULL;
double* eta = NULL;
double* m_j = NULL;

void kepler_step(int i){
	double M = G*(eta[i]);
	printf("\n %d %e\n",i,M);
	struct particle p1 = p_j[i];

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

	p_j[i].x = f*p1.x + g*p1.vx;
	p_j[i].y = f*p1.y + g*p1.vy;
	p_j[i].z = f*p1.z + g*p1.vz;

	p_j[i].vx = fd*p1.x + gd*p1.vx;
	p_j[i].vy = fd*p1.y + gd*p1.vy;
	p_j[i].vz = fd*p1.z + gd*p1.vz;

}


void integrator_part1(){
}


void integrator_to_jacobi(){
	p_j[0].x = 0.;
	p_j[0].y = 0.;
	p_j[0].z = 0.;
	for (int j=0;j<N;j++){
		p_j[0].x += 1./eta[N-1] * particles[j].m * particles[j].x;
		p_j[0].y += 1./eta[N-1] * particles[j].m * particles[j].y;
		p_j[0].z += 1./eta[N-1] * particles[j].m * particles[j].z;
	}
	for (int i=1;i<N;i++){
		p_j[i].x = particles[i].x;
		p_j[i].y = particles[i].y;
		p_j[i].z = particles[i].z;
		for (int j=0;j<i;j++){
			p_j[i].x -= 1./eta[i-1] * particles[j].m * particles[j].x;
			p_j[i].y -= 1./eta[i-1] * particles[j].m * particles[j].y;
			p_j[i].z -= 1./eta[i-1] * particles[j].m * particles[j].z;
		}
	}


	p_j[0].vx = 0.;
	p_j[0].vy = 0.;
	p_j[0].vz = 0.;
	for (int j=0;j<N;j++){
		p_j[0].vx += 1./eta[N-1] * particles[j].m * particles[j].vx;
		p_j[0].vy += 1./eta[N-1] * particles[j].m * particles[j].vy;
		p_j[0].vz += 1./eta[N-1] * particles[j].m * particles[j].vz;
	}
	for (int i=1;i<N;i++){
		p_j[i].vx = eta[i-1]/eta[i] * particles[i].m * particles[i].vx;
		p_j[i].vy = eta[i-1]/eta[i] * particles[i].m * particles[i].vy;
		p_j[i].vz = eta[i-1]/eta[i] * particles[i].m * particles[i].vz;
		for (int j=0;j<i;j++){
			p_j[i].vx -= particles[i].m/eta[i] * particles[j].m * particles[j].vx;
			p_j[i].vy -= particles[i].m/eta[i] * particles[j].m * particles[j].vy;
			p_j[i].vz -= particles[i].m/eta[i] * particles[j].m * particles[j].vz;
		}
		p_j[i].vx /= m_j[i];
		p_j[i].vy /= m_j[i];
		p_j[i].vz /= m_j[i];
	}
}

void integrator_to_heliocentric(){
	particles[0].x = p_j[0].x;
	particles[0].y = p_j[0].y;
	particles[0].z = p_j[0].z;
	for (int j=1;j<N;j++){
		particles[0].x -= particles[j].m/eta[j]*p_j[j].x;
		particles[0].y -= particles[j].m/eta[j]*p_j[j].y;
		particles[0].z -= particles[j].m/eta[j]*p_j[j].z;
	}
	for (int i=1;i<N;i++){
		particles[i].x = p_j[0].x + eta[i-1]/eta[i] * p_j[i].x;
		particles[i].y = p_j[0].y + eta[i-1]/eta[i] * p_j[i].y;
		particles[i].z = p_j[0].z + eta[i-1]/eta[i] * p_j[i].z;
		for (int j=i+1;j<N;j++){
			particles[i].x -= particles[j].m/eta[j] * p_j[j].x;
			particles[i].y -= particles[j].m/eta[j] * p_j[j].y;
			particles[i].z -= particles[j].m/eta[j] * p_j[j].z;
		}
	}
	particles[0].vx = particles[0].m/eta[N-1] * p_j[0].vx*m_j[0];
	particles[0].vy = particles[0].m/eta[N-1] * p_j[0].vy*m_j[0];
	particles[0].vz = particles[0].m/eta[N-1] * p_j[0].vz*m_j[0];
	for (int j=1;j<N;j++){
		particles[0].vx -= particles[0].m/eta[j-1]*p_j[j].vx*m_j[j];
		particles[0].vy -= particles[0].m/eta[j-1]*p_j[j].vy*m_j[j];
		particles[0].vz -= particles[0].m/eta[j-1]*p_j[j].vz*m_j[j];
	}
	particles[0].vx /= particles[0].m;
	particles[0].vy /= particles[0].m;
	particles[0].vz /= particles[0].m;

	for (int i=1;i<N;i++){
		particles[i].vx = particles[i].m/eta[N-1] * p_j[0].vx*m_j[0] + p_j[i].vx*m_j[i];
		particles[i].vy = particles[i].m/eta[N-1] * p_j[0].vy*m_j[0] + p_j[i].vy*m_j[i];
		particles[i].vz = particles[i].m/eta[N-1] * p_j[0].vz*m_j[0] + p_j[i].vz*m_j[i];
		for (int j=i+1;j<N;j++){
			particles[i].vx -= particles[i].m/eta[j-1]*p_j[j].vx*m_j[j];
			particles[i].vy -= particles[i].m/eta[j-1]*p_j[j].vy*m_j[j];
			particles[i].vz -= particles[i].m/eta[j-1]*p_j[j].vz*m_j[j];
		}
		particles[i].vx /= particles[i].m;
		particles[i].vy /= particles[i].m;
		particles[i].vz /= particles[i].m;
	}
}

void acceleration_jacobi(double _dt){
	for (int i=1;i<N;i++){
		// Eq 132
		double rj3 = pow(p_j[i].x*p_j[i].x+p_j[i].y*p_j[i].y+p_j[i].z*p_j[i].z,3./2.);
		p_j[i].vx += _dt * G*eta[i]/rj3 * p_j[i].x;
		p_j[i].vy += _dt * G*eta[i]/rj3 * p_j[i].y;
		p_j[i].vz += _dt * G*eta[i]/rj3 * p_j[i].z;
	}
}

void acceleration_heliocentric(double _dt){
	for (int i=1;i<N;i++){
		// Eq 132
		particles[i].ax = 0.;
		particles[i].ay = 0.;
		particles[i].az = 0.;
		for (int k=0;k<N;k++){
			if (k!=i){
				double rikx = particles[k].x-particles[i].x;
				double riky = particles[k].y-particles[i].y;
				double rikz = particles[k].z-particles[i].z;
				double rik3 = pow(rikx*rikx + riky*riky + rikz*rikz,3./2.);
				particles[i].ax += G*particles[k].m/rik3 * rikx;
				particles[i].ay += G*particles[k].m/rik3 * riky;
				particles[i].az += G*particles[k].m/rik3 * rikz;
			}
		}
		for (int j=0;j<i;j++){
			for (int k=0;k<N;k++){
				if (k!=j){
					double rjkx = particles[k].x-particles[j].x;
					double rjky = particles[k].y-particles[j].y;
					double rjkz = particles[k].z-particles[j].z;
					double rjk3 = pow(rjkx*rjkx + rjky*rjky + rjkz*rjkz,3./2.);
					particles[i].ax += 1./eta[i-1]*G*particles[j].m*particles[k].m/rjk3 * rjkx;
					particles[i].ay += 1./eta[i-1]*G*particles[j].m*particles[k].m/rjk3 * rjky;
					particles[i].az += 1./eta[i-1]*G*particles[j].m*particles[k].m/rjk3 * rjkz;
				}
			}
		}
		particles[i].vx += _dt * particles[i].ax;
		particles[i].vy += _dt * particles[i].ay;
		particles[i].vz += _dt * particles[i].az;
	}
}

void integrator_part2(){
	if (p_j==NULL){
		p_j = malloc(sizeof(struct particle)*N);
		eta = malloc(sizeof(double)*N);
		m_j = malloc(sizeof(double)*N);
		eta[0] = particles[0].m;
		for (int i=1;i<N;i++){
			eta[i] = eta[i-1] + particles[i].m;
		}
		m_j[0] = eta[N-1];
		for (int i=1;i<N;i++){
			m_j[i] = eta[i-1]/eta[i]*particles[i].m;
		}
	}




	acceleration_heliocentric(dt/2.);
	integrator_to_jacobi();
	acceleration_jacobi(dt/2.);

	for (int i=1;i<N;i++){
		kepler_step(i);
	}
	p_j[0].x += dt*p_j[0].vx;
	p_j[0].y += dt*p_j[0].vy;
	p_j[0].z += dt*p_j[0].vz;

	acceleration_jacobi(dt/2.);
	integrator_to_heliocentric();
	acceleration_heliocentric(dt/2.);

	t+=dt;
}
	

