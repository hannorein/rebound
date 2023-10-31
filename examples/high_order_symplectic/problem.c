/**
 * High Order Symplectic Integrators
 *
 * This example uses a high order symplectic integrators
 * WHCKL and SABA(10,6,4) to integrate all planets of the Solar System. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double ss_pos[10][3] = 
{
    {3.256101656448802E-03  , -1.951205394420489E-04 , -1.478264728548705E-04},  
    {-1.927589645545195E-01 , 2.588788361485397E-01  , 3.900432597062033E-02 }, 
    {-5.976537074581466E-01 , 3.918678996109574E-01  , 3.990356741282203E-02 }, 
    {-7.986189029000561E-01 , -6.086873314992410E-01 , -1.250824315650566E-04}, 
    {7.897942807177173E-01  , 1.266671734964037E+00  , 7.092292179885432E-03 }, 
    {-4.314503046344270E+00 , 3.168094294126697E+00  , 8.331048545353310E-02 }, 
    {-4.882304833383455E+00 , -8.689263067189865E+00 , 3.453930436208210E-01 }, 
    {1.917757033372740E+01  , 5.671738750949031E+00  , -2.273858614425555E-01},  
    {2.767031517959636E+01  , -1.150331645280942E+01 , -4.008018419157927E-01},  
    {7.765250227278298E+00  , -3.190996242617413E+01 , 1.168394015703735E+00 }, 

};
double ss_vel[10][3] = 
{
    {3.039963463108432E-06 ,  6.030576499910942E-06 ,  -7.992931269075703E-08}, 
    {-2.811550184725887E-02,  -1.586532995282261E-02,  1.282829413699522E-03 }, 
    {-1.113090630745269E-02,  -1.703310700277280E-02,  4.089082927733997E-04 },
    {1.012305635253317E-02 ,  -1.376389620972473E-02,  3.482505080431706E-07 }, 
    {-1.135279609707971E-02,  8.579013475676980E-03 ,  4.582774369441005E-04 }, 
    {-4.555986691913995E-03,  -5.727124269621595E-03,  1.257262404884127E-04 }, 
    {4.559352462922572E-03 ,  -2.748632232963112E-03,  -1.337915989241807E-04}, 
    {-1.144087185031310E-03,  3.588282323722787E-03 ,  2.829006644043203E-05 }, 
    {1.183702780101068E-03 ,  2.917115980784960E-03 ,  -8.714411604869349E-05}, 
    {3.112825364672655E-03 ,  1.004673400082409E-04 ,  -9.111652976208292E-04},
};

double ss_mass[10] =
{
    1.988544e30,
    3.302e23,
    48.685e23,
    6.0477246e24,
    6.4185e23,
    1898.13e24,
    5.68319e26,
    86.8103e24,
    102.41e24,
    1.4639248e+22,
};

struct reb_simulation* create_sim(){
    // Setup constants
    struct reb_simulation* r = reb_simulation_create();
    r->dt             = 4;          // in days
    r->G              = 1.0e-34;    // in AU^3 / kg / day^2.
    
    // Initial conditions (from NASA Horizons)
    for (int i=0;i<10;i++){
        struct reb_particle p = {0};
        p.x  = ss_pos[i][0];         p.y  = ss_pos[i][1];         p.z  = ss_pos[i][2];
        p.vx = ss_vel[i][0];         p.vy = ss_vel[i][1];         p.vz = ss_vel[i][2];
        p.m  = ss_mass[i];
        reb_simulation_add(r, p); 
    }

    return r;
}

int main(int argc, char* argv[]){
    double tmax       = 1e5;              // 1e5 days ~ 273 years

    // Run the simulation with the WHCKL method.
    {    
        struct reb_simulation* r = create_sim();
        r->integrator           = REB_INTEGRATOR_WHFAST;
        r->ri_whfast.safe_mode  = 0;        // Turn off safe mode (Need to call reb_simulation_synchronize() before outputs).
        r->ri_whfast.corrector  = 17;       // 17th order symplectic corrector
        r->ri_whfast.kernel     = REB_WHFAST_KERNEL_LAZY;   // Using the lazy implementers method which supports additional forces
        double e_init = reb_simulation_energy(r);
        reb_simulation_integrate(r, tmax);
        double e = reb_simulation_energy(r);
        printf("Relative energy error WHCKL: %e\n", fabs((e_init-e)/e_init));
    }
    
    // Run the same simulation with the SABA(10,6,4) method.
    // Note that this method has 8 force evaluations per timestep and is therefore 
    // quite a bit slower for a fixed timestep.
    {    
        struct reb_simulation* r = create_sim();
        r->integrator           = REB_INTEGRATOR_SABA; 
        r->ri_saba.type  = REB_SABA_10_6_4;    // Chooses the type of SABA integrator. 
        r->ri_saba.safe_mode  = 0;        // Turn off safe mode. 
        double e_init = reb_simulation_energy(r);
        reb_simulation_integrate(r, tmax);
        double e = reb_simulation_energy(r);
        printf("Relative energy error SABA(10,6,4):    %e\n", fabs((e_init-e)/e_init));
    }
    // Run the same simulation with the standard WH method.
    {    
        struct reb_simulation* r = create_sim();
        r->integrator           = REB_INTEGRATOR_WHFAST; // All WHFast settings default to the standard WH method
        r->ri_whfast.safe_mode  = 0;        // Turn off safe mode. 
        double e_init = reb_simulation_energy(r);
        reb_simulation_integrate(r, tmax);
        double e = reb_simulation_energy(r);
        printf("Relative energy error WH:    %e\n", fabs((e_init-e)/e_init));
    }
}

