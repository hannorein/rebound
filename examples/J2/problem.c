/**
 * J2 precession
 *
 * This example presents an implementation of the J2 gravitational moment.
 * The equation of motions are integrated with the 15th order IAS15
 * integrator. The parameters in this example have been chosen to
 * represent those of Saturn, but one can easily change them or even
 * include higher order terms in the multipole expansion.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

const double J2planet           = 16298e-6;         // J2 of Saturn (Murray and Dermott p 531)
const double Mplanet            = 0.00028588598;    // mass of Saturn in solar masses
const double Rplanet            = 0.00038925688;    // radius of Saturn in AU
const double ObliquityPlanet    = 0.;               // obliquity of the planet

const double tmax           = 1e3;          // Maximum integration time

void heartbeat(struct reb_simulation* r);
void force_J2(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator           = REB_INTEGRATOR_IAS15;
    r->dt               = 1e-6;         // initial timestep
    r->N_active         = 2;            // only the star and the planet are massive.

    // Planet
    struct reb_particle planet = {0};
    planet.m  = Mplanet;
    reb_simulation_add(r, planet);

    struct reb_particle p = {0};          // test particle
    double a = Rplanet*3.;                // small distance from planet (makes J2 important)
    double e = 0.1;
    double v = sqrt((1.+e)/(1.-e)*r->G*planet.m/a); // setup eccentric orbit (ignores J2)
    p.x  = (1.-e)*a;
    p.vy = v;
    p.x += planet.x;    p.y += planet.y;    p.z += planet.z;
    p.vx += planet.vx;  p.vy += planet.vy;  p.vz += planet.vz;
    reb_simulation_add(r, p);

    reb_simulation_move_to_com(r);

    remove("a.txt");                      // delete previous output

    // Setup callback functions
    r->heartbeat        = heartbeat;
    r->additional_forces    = force_J2;

    reb_simulation_integrate(r, tmax);

    reb_simulation_free(r);
}

void force_J2(struct reb_simulation* r){
    if (J2planet==0) return;
    // Star
    const struct reb_particle planet = r->particles[0];     // cache
    const int N = r->N;
#pragma omp parallel for
    for (int i=1;i<N;i++){
        const struct reb_particle p = r->particles[i];      // cache
        const double sprx = p.x-planet.x;
        const double spry = p.y-planet.y;
        const double sprz = p.z-planet.z;
        const double prx  = sprx*cos(-ObliquityPlanet) + sprz*sin(-ObliquityPlanet);
        const double pry  = spry;
        const double prz  =-sprx*sin(-ObliquityPlanet) + sprz*cos(-ObliquityPlanet);
        const double pr2  = prx*prx + pry*pry + prz*prz;        // distance^2 relative to planet
        const double fac  = 3.*r->G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);

        const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
        const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
        const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);

        r->particles[i].ax += pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet);
        r->particles[i].ay += pay;
        r->particles[i].az +=-pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet);
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r, 4000.*r->dt)){                // output something to screen
        reb_simulation_output_timing(r, tmax);
    }
    if(reb_simulation_output_check(r,M_PI*2.*0.01)){                // output some orbital parameters to file
        FILE* f = fopen("a.txt","ab");
        const struct reb_particle planet = r->particles[0];
        const int N = r->N;
        for (int i=1;i<N;i++){
            struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[i],planet);
            fprintf(f,"%.15e\t%.15e\t%.15e\t%.15e\n",r->t,o.a,o.e,o.omega);
        }
        fclose(f);
    }
}

