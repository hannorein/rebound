/**
 * The chaos indicator MEGNO
 * 
 * This example uses the IAS15 or WHFAST integrator
 * to calculate the MEGNO of a two planet system.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

const double ss_pos[3][3] = 
{
    {-4.06428567034226e-3,    -6.08813756435987e-3,    -1.66162304225834e-6    }, // Sun
    {+3.40546614227466e+0,    +3.62978190075864e+0,    +3.42386261766577e-2    }, // Jupiter
    {+6.60801554403466e+0,    +6.38084674585064e+0,    -1.36145963724542e-1    }, // Saturn
};
const double ss_vel[3][3] = 
{
    {+6.69048890636161e-6,    -6.33922479583593e-6,    -3.13202145590767e-9    }, // Sun
    {-5.59797969310664e-3,    +5.51815399480116e-3,    -2.66711392865591e-6    }, // Jupiter
    {-4.17354020307064e-3,    +3.99723751748116e-3,    +1.67206320571441e-5    }, // Saturn
};

const double ss_mass[3] =
{
    1.00000597682,     // Sun + inner planets
    1./1047000.355,    // Jupiter
    1./3501000.6,    // Saturn
};

double tmax = 1e9;

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    // Setup constants
    r->dt             = 10;            // initial timestep (in days)
    //r->integrator    = IAS15;
    r->integrator     = REB_INTEGRATOR_WHFAST;
    const double k    = 0.01720209895; // Gaussian constant 
    r->G              = k*k;           // These are the same units that mercury6 uses

    // Initial conditions
    for (int i=0;i<3;i++){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];
        reb_simulation_add(r, p); 
    }
    reb_simulation_move_to_com(r);
    // Add megno particles 
    reb_simulation_init_megno(r);  // N = 6 after this function call. 
    // The first half of particles are real particles, the second half are particles following the variational equations.
    
    // Set callback for outputs.
    r->heartbeat = heartbeat;

    reb_simulation_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 100.*362.)){
        // Output the time and the MEGNO to the screen and a file every 100 years.
        FILE* f = fopen("Y.txt","a+b");
        fprintf(f,"        %.20e     %.20e\n",r->t, reb_simulation_megno(r));
        printf(" t= %.2e   MEGNO = %.2e\n",r->t, reb_simulation_megno(r));
        fclose(f);
    }
}

void problem_finish(){
}
