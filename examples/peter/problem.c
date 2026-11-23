#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include "chariklo.h"
#define N_p 24000 // ring particle number

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);
void heartbeat(struct reb_simulation* const r);
void force_J2_C22(struct reb_simulation* r);

const int tmax = 10000; // time during which impact velocities are measured (measured in periods)
const double dt_frac = 1./1200.;
struct reb_timeval tim_start;

int main(int argc, char* argv[]) {
    // Get the number of processors
    int np = omp_get_num_procs();
    printf("Number of processors: %d\n", np);
    // Set the number of OpenMP threads to be the number of processors    
    //omp_set_num_threads(np);
    #pragma omp parallel
    {
        #pragma omp single
        printf("Inside parallel region: %d\n", omp_get_num_threads());
    }
	struct reb_simulation* r = reb_simulation_create();
    reb_simulation_start_server(r, 1234);
	// Setup constants
	reb_simulation_set_integrator(r, "leapfrog");
    r->gravity      = REB_GRAVITY_TREE;
    r->collision    = REB_COLLISION_LINETREE;
    r->boundary     = REB_BOUNDARY_OPEN;
    r->opening_angle2 = .5;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;

    double boxsize = 3.5e7;
    r->root_size = boxsize;
    r->N_root_x = 2;
    r->N_root_y = 2;
	
    r->G                 = 6.67428e-11;         // N / kg^2 m^2
    double a_cor = cbrt(r->G*Mplanet/pow(OMEGA_H, 2)); // corotation radius
    double a0 = cbrt(r->G*Mplanet/pow(OMEGA_H/3., 2)); // radius corresponding to 1/3 SOR
    double OMEGA = sqrt(r->G*Mplanet/pow(a0, 3));
	r->dt                = 2*M_PI/OMEGA*dt_frac;  // s
	r->heartbeat         = heartbeat;           // function pointer for heartbeat
   	double particle_density = 300; // kg/m^3
    double particle_radius = 0.005*a_cor;
    struct reb_particle planet = {0};
    planet.m  = Mplanet;
    planet.r = Rplanet;
    reb_simulation_add(r, planet);
    
    double a, m, true_anomaly;
    for (int i = 0; i < N_p; i++) {
        a = reb_random_uniform(r, 1.55*a_cor, 2.8*a_cor);
        true_anomaly = reb_random_uniform(r, 0., 2*M_PI);
        m = 4./3.*particle_radius*particle_radius*particle_radius*M_PI*particle_density;
        reb_simulation_add_fmt(r, "r a f m", particle_radius, a, true_anomaly, m);
    }

    r->additional_forces = force_J2_C22;
    const double period_b = 2*M_PI/OMEGA_H; // rotational period of central body
		
    gettimeofday(&tim_start, NULL);
    reb_simulation_integrate(r, tmax*period_b);

    reb_simulation_free(r);
}

/// @brief Zonal and tesseral contributions are computed in the body-fixed, rotating frame and then converted back
/// into the inertial frame.
/// @param r 
void force_J2_C22(struct reb_simulation* r){
    const struct reb_particle planet = r->particles[0];
    double j2 = J2planet, c22 = C22planet;

    if (planet.m == 0. || (j2 == 0. && c22 == 0.)) return;

    const double angle = c22 == 0. ? 0. : OMEGA_H*r->t;
    const double c = cos(angle);
    const double s = sin(angle);
    const double K = -3.*r->G*planet.m*Rplanet*Rplanet;
    const int N = r->N;
    double Fx = 0., Fy = 0., Fz = 0.;

    #pragma omp parallel for reduction(+:Fx,Fy,Fz)
    for (int i = 1; i < N; i++){
        struct reb_particle* p = &r->particles[i];
        const double dx = p->x - planet.x;
        const double dy = p->y - planet.y;
        const double z = p->z - planet.z;
        const double x = c*dx + s*dy;
        const double y = -s*dx + c*dy;
        const double x2 = x*x, y2 = y*y, z2 = z*z;
        const double invr2 = 1./(x2 + y2 + z2);
        const double fac = K*invr2*invr2*invr2*sqrt(invr2);
        const double bx = fac*x*(
            0.5*j2*(x2 + y2 - 4.*z2)
            + c22*(2.*z2 + 7.*y2 - 3.*x2));
        const double by = fac*y*(
            0.5*j2*(x2 + y2 - 4.*z2)
            + c22*(3.*y2 - 2.*z2 - 7.*x2));
        const double bz = fac*z*(
            0.5*j2*(3.*(x2 + y2) - 2.*z2)
            + 5.*c22*(y2 - x2));
        const double ax = c*bx - s*by;
        const double ay = s*bx + c*by;
        const double az = bz;
        p->ax += ax;
        p->ay += ay;
        p->az += az;
        Fx += p->m*ax;
        Fy += p->m*ay;
        Fz += p->m*az;
    }
    r->particles[0].ax -= Fx/planet.m;
    r->particles[0].ay -= Fy/planet.m;
    r->particles[0].az -= Fz/planet.m;
}

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void heartbeat(struct reb_simulation* const r){
	if (reb_simulation_output_check(r, 1e-1*2.*M_PI/OMEGA_H)){
		struct reb_timeval tim;
        gettimeofday(&tim, NULL);
        double temp = tim.tv_sec+(tim.tv_usec/1000000.0) - tim_start.tv_sec+(tim_start.tv_usec/1000000.0);
        tim_start = tim;
        printf("N_tot= %- 9ld  ",r->N);
        printf("t= %- 9f [orb]  ",r->t*OMEGA_H/2./M_PI);
        printf("dt= %- 9f  ",r->dt);
        printf("cpu= %- 9f [s]  \n",temp);
	}
}


