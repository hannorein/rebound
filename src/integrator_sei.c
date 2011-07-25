#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "gravity.h"
#include "main.h"
#include "boundaries.h"

double OMEGA = 1.; // Epicyclic frequency 
double OMEGAZ = -1.; // Epicyclic frequency in vertical direction

void operator_H0(double dt, struct particle* p);
void operator_phi(double dt, struct particle* p);
// Cache sin() tan() values.
double lastdt=0;
double sindt, tandt;
double sindtz, tandtz;

// This function is the SEI integrator.
// It is symplectic, second order accurate and time-reversible.
// I.e. pretty cool.
void integrator_part1(){
	if (lastdt!=dt){
		// Only calculate sin() and tan() if timestep changed
		if (OMEGAZ==-1){
			OMEGAZ=OMEGA;
		}
		sindt = sin(OMEGA*(-dt));
		tandt = tan(OMEGA*(-dt/2.));
		sindtz = sin(OMEGAZ*(-dt));
		tandtz = tan(OMEGAZ*(-dt/2.));
		lastdt = dt;
	}

#pragma omp parallel for
	for (int i=0;i<N;i++){
		operator_H0(dt/2.,&(particles[i]));
	}
	t+=dt/2.;
}

void integrator_part2(){
#pragma omp parallel for
	for (int i=0;i<N;i++){
		operator_phi(dt,&(particles[i]));
		operator_H0(dt/2.,&(particles[i]));
	}
	t+=dt/2.;
}

// This function evolves a particle under 
// Hamiltonian H0 exactly up to machine precission.
void operator_H0(double dt, struct particle* p){
		
	// Integrate vertical motion
	const double zx = p->z * OMEGAZ;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - tandtz*zy;			
	const double zyt =  sindtz*zt1 + zy;
	const double zxt =  zt1 - tandtz*zyt;	
	p->z  = zxt/OMEGAZ;
	p->vz = zyt;

	// Integrate motion in xy directions
	const double aO = 2.*p->vy + 4.*p->x*OMEGA;	// Center of epicyclic motion
	const double bO = p->y*OMEGA - 2.*p->vx;	

	const double ys = (p->y*OMEGA-bO)/2.; 		// Epicycle vector
	const double xs = (p->x*OMEGA-aO); 
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double xst1 =  xs - tandt*ys;			
	const double yst  =  sindt*xst1 + ys;
	const double xst  =  xst1 - tandt*yst;	

	p->x  = (xst+aO)    /OMEGA;			
	p->y  = (yst*2.+bO) /OMEGA - 3./2.*aO*dt;	
	p->vx = yst;
	p->vy = -xst*2. -3./2.*aO;
}

// This function evolves a particle under the operator
// Phi exactly up to machine precission.
void operator_phi(double dt, struct particle* p){
	// The force used here is for test cases 2 and 3 
	// in Rein & Tremaine 2011. 
	p->vx += p->ax * dt;
	p->vy += p->ay * dt;
	p->vz += p->az * dt;
}

