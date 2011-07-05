#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "gravity.h"
#include "main.h"

const double OMEGA = 1.; // Orbital velocity

void operator_H0(double dt, struct particle* p);
void operator_phi(double dt, struct particle* p);
// Cache sin() tan() values.
double lastdt=0;
double sindt, tandt;

// This function is the SEI integrator.
// It is symplectic, second order accurate and time-reversible.
// I.e. pretty cool.
void integrate_particles(){
	for (int i=0;i<N;i++){
		operator_H0(dt/2.,&(particles[i]));
	}
	calculate_forces();
	for (int i=0;i<N;i++){
		operator_phi(dt,&(particles[i]));
		operator_H0(dt/2.,&(particles[i]));
	}
}

// This function evolves a particle under 
// Hamiltonian H0 exactly up to machine precission.
void operator_H0(double dt, struct particle* p){
	if (lastdt!=dt){
		// Only calculate sin() and tan() if timestep changed
		sindt = sin(OMEGA*(-dt));
		tandt = tan(OMEGA*(-dt/2.));
		lastdt = dt;
	}
		
	// Integrate vertical motion
	const double zx = p->z * OMEGA;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - tandt*zy;			
	const double zyt =  sindt*zt1 + zy;
	const double zxt =  zt1 - tandt*zyt;	
	p->z  = zxt/OMEGA;
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

