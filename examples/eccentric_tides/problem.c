/*
*
*
*  Tidal decay of eccentric orbits in two body problem
*
*  The tidal buldge rises as per eccentric orbit, which 
*  induces a const. lag as due to the asyn. rotation of   
*  the planet and its parent star. 
*
*  As the total angular momentum conserves, while energy 
*  dissipates continuesly. For a 3-day (~ 0.04 au) Hot
*  Jupiter - Sun-like star system, Keplerian Orbital AM 
*  ~ 1.7e+49 [cgs units]; The rotational momentum of the 
*  central star (Prot ~ 20 days ) is around 1.4e+49 - they
*  are comparable! However, the distortion on the planet 
*  caused by the star is more significant than the couter 
*  part raised on the star.
*
*
*  by Meldonization. 2016
*
*
                   _ooOoo_
                  o8888888o
                  88" . "88
                  (| -_- |)
                  O\  =  /O
               ____/`---'\____
             .'  \\|     |//  `.
            /  \\|||  :  |||//  \
           /  _||||| -:- |||||-  \
           |   | \\\  -  /// |   |
           | \_|  ''\---/''  |   |
           \  .-\__  `-`  ___/-. /
         ___`. .'  /--.--\  `. . __
      ."" '<  `.___\_<|>_/___.'  >'"".
     | | :  `- \`.;`\ _ /`;.`/ - ` : | |
     \  \ `-.   \_ __\ /__ _/   .-` /  /
======`-.____`-.___\_____/___.-`____.-'======
                   `=---='
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          God Bless       Zero Bug

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void tidal_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);
double tmax;
double e_init; // initial energy

int main(int argc, char* argv[]) { 
	struct reb_simulation* r = reb_create_simulation();
	tmax	= 1e8 * 2.* M_PI;
	
	// Setup constants
	r->dt 			= 1e-3;		// initial timestep.
	r->G			= 1;		// Gravitational constant
	r->integrator	= REB_INTEGRATOR_IAS15;

	// Setup callback function for velocity dependent forces.
	r->additional_forces 	= tidal_forces;
	r->force_is_velocity_dependent = 1;
	// Setup callback function for outputs.
	r->heartbeat	= heartbeat;
	r->usleep		= 10000;		
	double mass_scale	= 1.;		
	double size_scale	= 1.;		
	double rads_scale	= 5.0e-3;	// Solar radius in au.	
	double e_jup	 	= 0.985;	
	
	struct reb_particle star = {0}; 
	star.m = mass_scale;
	star.r = rads_scale;
	reb_add(r, star); 
	
	struct reb_particle p; 
	p.m  = mass_scale * 9.54265748e-4;	// Jupiter mass in M_sun
	p.r  = rads_scale * 0.1;  // Jupiter radius in R_sun
	p.x  = size_scale * (1. + e_jup); 
	p.vy = sqrt((1. - e_jup) / (1. + e_jup) * mass_scale / size_scale);
	reb_add(r, p); 

	e_init = reb_tools_energy(r);
	reb_move_to_com(r);
	system("rm -v orbits.txt"); // delete previous output file
	
	// Do the integration
	reb_integrate(r, tmax);

}

static double acos2(double num, double denom, double disambiguator){
	double val;
	double cosine = num/denom;
	if(cosine > -1. && cosine < 1.){
		val = acos(cosine);
		if(disambiguator < 0.){
			val = - val;
		}
	}
	else{
		val = (cosine <= -1.) ? M_PI : 0.;
	}
	return val;
}

const double eccsmooth=1.E-5;
/* tiny eccentricity causes repeat calculating true anomaly problem */
int ecctrigger = 0;

void tidal_forces(struct reb_simulation* const r){	
	struct reb_particle* const particles = r->particles;
	struct reb_particle com = particles[0];
	const double Qp = 1.0e5; // tidal deformation parameter
	const double k2p = 5.2e-1; // love number
	const double G = r->G;
	const int N = r->N;
	double coef, coefr, taolag, rdot, f;
	
	for (int i=1;i<N;i++){
		
		if (ecctrigger) {
			continue;
		}
		
		struct reb_particle* p = &(particles[i]);
		const double dx = p->x-com.x;
		const double dy = p->y-com.y;
		const double dz = p->z-com.z;
		const double d = sqrt(dx*dx + dy*dy + dz*dz);
		
		//if (d >= 0.5 ) continue; // tidal force effective at close distances
		
		const double dvx = p->vx-com.vx;
		const double dvy = p->vy-com.vy;
		const double dvz = p->vz-com.vz;
		//const double hx = dy*dvz - dz*dvy; 
		//const double hy = dz*dvx - dx*dvz;
		//const double hz = dx*dvy - dy*dvx;
		const double mu = G*(com.m + p->m);
		//const double h = sqrt(hx*hx + hy*hy + hz*hz);
		const double v = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
		const double vr = (dx*dvx + dy*dvy + dz*dvz)/d;
		const double ex = 1./mu*((v*v-mu/d)*dx - d*vr*dvx);
		const double ey = 1./mu*((v*v-mu/d)*dy - d*vr*dvy);
		const double ez = 1./mu*((v*v-mu/d)*dz - d*vr*dvz);
		const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
		
		if (e <= eccsmooth) { 
			// small eccentricity would slow down the intergration.
			ecctrigger = 1;
		} else if (e >= eccsmooth * 100) {
			ecctrigger = 0;
		}
		
		const double a = -mu/( v*v - 2.*mu/d );	
		const double n = sqrt(mu / pow(a, 3.0));
		const double q = a * (1. - e * e);
				
		f = acos2(q/d-1., e, vr);
		rdot = n * a * e / sqrt(1.0e0 - e * e) * sin(f);
		taolag = M_PI / Qp / n;
		coef = mu * com.m / p->m * k2p * pow(p->r, 5.0e0) / pow(d, 8.0e0);
		coefr = 9.0e0 * rdot / d * taolag;
	
		p->ax -= coef * coefr * dx;
		p->ay -= coef * coefr * dy;
		p->az -= coef * coefr * dz;
		
		reb_move_to_com(r); // avoid com drifting
	}		
		
}


void heartbeat(struct reb_simulation* const r){
		reb_output_timing(r, tmax);
		
		if(reb_output_check(r, M_PI*2000.)){
			//reb_integrator_synchronize(r);
			double en = reb_tools_energy(r);
			const struct reb_particle* p = r->particles;
			const double G = r->G;
			const int N = r->N;
			
			for (int i=1;i<N;i++){
				
				const double dx = p[i].x-p[0].x;
				const double dy = p[i].y-p[0].y;
				const double dz = p[i].z-p[0].z;
				const double d = sqrt(dx*dx + dy*dy + dz*dz);
				const double dvx = p[i].vx-p[0].vx;
				const double dvy = p[i].vy-p[0].vy;
				const double dvz = p[i].vz-p[0].vz;
				const double mu = G*(p[0].m + p[i].m);
				const double v = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
				const double vr = (dx*dvx + dy*dvy + dz*dvz)/d;
				const double ex = 1./mu*((v*v-mu/d)*dx - d*vr*dvx);
				const double ey = 1./mu*((v*v-mu/d)*dy - d*vr*dvy);
				const double ez = 1./mu*((v*v-mu/d)*dz - d*vr*dvz);
				const double e = sqrt( ex*ex + ey*ey + ez*ez );	// eccentricity
				const double a = -mu/( v*v - 2.*mu/d );	
				const double q = a * (1. - e);
				//struct reb_orbit o = reb_tools_particle_to_orbit(r->G, p[1], p[0]);
			
				FILE* f = fopen("orbits.txt","a");
				fprintf(f,"%e, %e, %e, %e, %e\n", r->t, a, q, e, fabs((en-e_init)/e_init));
				fclose(f);
			}
		}
}
