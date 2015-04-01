#include "librebound.h"
#include "tools.h"
#include "particle.h"
#include <stdlib.h>
#include <math.h>

int DISK_FORCES = 0;

struct disk{
	double gamma;
	double Rc;
	double m;
	double alpha;
	double podot; // pericenter precession at r = Rc
};

struct disk dsk; // structure to hold disk parameters for pericenter precession

// pointers for damping timescales
double *tau_a = NULL;
double *tau_e = NULL;
double *tau_i = NULL;

double e_damping_p; // p parameter from Goldreich & Schlichting 2014 for how e-damping
// contributes to a-damping at order e^2
// p = 3 : e-damping at const angular momentum.  p = 0 : no contribution to a-damping

void disk_forces(){
	struct particle com = particles[0]; // calculate add. forces w.r.t. center of mass
	for(int i=1;i<N;i++){
		struct particle* p = &(particles[i]);
		const double dvx = p->vx - com.vx;
		const double dvy = p->vy - com.vy;
		const double dvz = p->vz - com.vz;

		if (tau_a[i] != 0.){
			p->ax -=  dvx/(2.*tau_a[i]);
			p->ay -=  dvy/(2.*tau_a[i]);
			p->az -=  dvz/(2.*tau_a[i]);
		}

		if (tau_e[i] != 0. || tau_i[i]!= 0. || dsk.m != 0.){ 	// need h and e vectors for both types
			const double mu = G*(com.m + p->m);
			const double dx = p->x-com.x;
			const double dy = p->y-com.y;
			const double dz = p->z-com.z;

			const double hx = dy*dvz - dz*dvy;
			const double hy = dz*dvx - dx*dvz;
			const double hz = dx*dvy - dy*dvx;
			const double h = sqrt ( hx*hx + hy*hy + hz*hz );
			const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
			const double r = sqrt ( dx*dx + dy*dy + dz*dz );
			const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
			const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
			const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
			const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
			const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity

			if (tau_e[i] != 0.){	// Eccentricity damping
				const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
				const double prefac1 = 1./tau_e[i]/1.5*(1.+e_damping_p/2.*e*e);
				const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))/tau_e[i]/1.5;

				p->ax += -2/tau_e[i]*vr*dx/r;
				p->ay += -2/tau_e[i]*vr*dy/r;
				p->az += -2/tau_e[i]*vr*dz/r;
				/*p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
				p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
				p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;*/
			}
			if (tau_i[i]!=0){		// Inclination damping
				p->az += -2.*dvz/tau_i[i];
				const double prefac = (hx*hx + hy*hy)/h/h/tau_i[i];
				p->ax += prefac*dvx;
				p->ay += prefac*dvy;
				p->az += prefac*dvz;
			}
			if (dsk.m != 0.) {
				double a_over_r = -dsk.alpha*pow(dsk.Rc/r,dsk.gamma)/r + G*dsk.m/r/r/r; 	// radial disk force after removing piece from adding the disk into the sun
				p->ax += a_over_r*dx;									// rhat has components x/r xhat + y/r yhat + z/r zhat
				p->ay += a_over_r*dy;
				p->az += a_over_r*dz;

				particles[0].ax -= p->m/particles[0].m*a_over_r*dx;		// add back reactions onto the star (if forces are equal, accelerations differ by -mass ratio)
				particles[0].ay -= p->m/particles[0].m*a_over_r*dy;
				particles[0].az -= p->m/particles[0].m*a_over_r*dz;
			}
		}
		com = tools_get_center_of_mass(com,particles[i]);
	}
	tools_move_to_center_of_momentum();
}

void add_migration(double *_tau_a){
	if(!DISK_FORCES) { // have to set other taus to 0 so can access when checking in force calc whether we should calculate them
		tau_a = calloc(sizeof(double),N);
		tau_e = calloc(sizeof(double),N);
		tau_i = calloc(sizeof(double),N);
	}

	for(int i=0; i<N; ++i){
		tau_a[i] = _tau_a[i];
	}
}

void set_e_damping(double *_tau_e){
	for(int i=0; i<N; ++i){
		tau_e[i] = _tau_e[i];
	}
}
void add_e_damping(double *_tau_e){
	if(!DISK_FORCES) {
		tau_a = calloc(sizeof(double),N);
		tau_e = calloc(sizeof(double),N);
		tau_i = calloc(sizeof(double),N);
	}

	set_e_damping(_tau_e);
}

void add_i_damping(double *_tau_i){
	if(!DISK_FORCES) {
		tau_a = calloc(sizeof(double),N);
		tau_e = calloc(sizeof(double),N);
		tau_i = calloc(sizeof(double),N);
	}

	for(int i=0; i<N; ++i){
		tau_i[i] = _tau_i[i];
	}
}

void add_peri_precession(double gam, double Rc, double podot){
	if(!DISK_FORCES) {
		tau_a = calloc(sizeof(double),N);
		tau_e = calloc(sizeof(double),N);
		tau_i = calloc(sizeof(double),N);
	}

	dsk.Rc = Rc;
	dsk.gamma = gam;
	dsk.podot = podot; // as a fraction of the mean motion

	dsk.alpha = G*particles[0].m/Rc/Rc*podot/(1.-gam/2.);
	dsk.m = 3.65557*particles[0].m*podot;
	particles[0].m += dsk.m;

}

void migration_reset(){
	free(tau_a);
	tau_a = NULL;
	free(tau_e);
	tau_e = NULL;
	free(tau_i);
	tau_i = NULL;
}


