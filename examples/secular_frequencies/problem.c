/**
 * Secular Frequencies
 *
 * This example integrates the outer Solar System and then performs a 
 * frequency analysis using the Frequency Modified Fourier Transform
 * to determine the secular frequencies (g-modes).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double ss_pos[6][3] =
    {
     {-4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6}, // Sun
     {+3.40546614227466e+0, +3.62978190075864e+0, +3.42386261766577e-2}, // Jupiter
     {+6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1}, // Saturn
     {+1.11636331405597e+1, +1.60373479057256e+1, +3.61783279369958e-1}, // Uranus
     {-3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1}, // Neptune
     {-2.13858977531573e+1, +3.20719104739886e+1, +2.49245689556096e+0}  // Pluto
};
double ss_vel[6][3] =
    {
     {+6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9}, // Sun
     {-5.59797969310664e-3, +5.51815399480116e-3, -2.66711392865591e-6}, // Jupiter
     {-4.17354020307064e-3, +3.99723751748116e-3, +1.67206320571441e-5}, // Saturn
     {-3.25884806151064e-3, +2.06438412905916e-3, -2.17699042180559e-5}, // Uranus
     {-2.17471785045538e-4, -3.11361111025884e-3, +3.58344705491441e-5}, // Neptune
     {-1.76936577252484e-3, -2.06720938381724e-3, +6.58091931493844e-4}  // Pluto
};

double ss_mass[6] =
    {
     1.00000597682, // Sun + inner planets
     1. / 1047.355, // Jupiter
     1. / 3501.6,   // Saturn
     1. / 22869.,   // Uranus
     1. / 19314.,   // Neptune
     0.0  // Pluto
};

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    const double k = 0.01720209895; // Gaussian constant
    r->dt = 120;                    // Timestep is 120 days.
    r->G = k * k;                   // These are the same units as used by the mercury6 code.
    r->integrator = REB_INTEGRATOR_WHFAST;

    // Initial conditions
    for (int i = 0; i < 6; i++) {
        struct reb_particle p = {0};
        p.x = ss_pos[i][0];
        p.y = ss_pos[i][1];
        p.z = ss_pos[i][2];
        p.vx = ss_vel[i][0];
        p.vy = ss_vel[i][1];
        p.vz = ss_vel[i][2];
        p.m = ss_mass[i];
        reb_simulation_add(r, p);
    }

    reb_simulation_move_to_com(r);
    
    int Nsamples = 2048; // Number of samples. Must be a power of two.
                         // Choose a larger number for better accuracy, e.g. 32768.
    double* inp = malloc(sizeof(double)*2*Nsamples);
    // Start integration
    for (int i=0; i<Nsamples; i++){
        // Integrate for 1000 steps (120000 days)
        reb_simulation_steps(r, 1000); 
        // Calculate orbital elements of Jupiter
        struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
        // Store complex eccentricity in array
        inp[i*2+0] = o.e*cos(o.pomega);
        inp[i*2+1] = o.e*sin(o.pomega);
    }

    // Perform frequency analysis
    int nfreq = 5;
    double datasep = 120000.0/365.25*2.0*M_PI;  // sampling interval in units of year/2pi
    double minfreq = 60.0/1296000.0*datasep;    // min/max frequenxy 60"/year 
    double* out = malloc(sizeof(double)*3*nfreq);
    // The next command performs the actual Frequency Modified Fourier Transform (FMFT). 
    // Other options are MFT (faster) and FMFT2 (more accurate). 
    // See Sidlichovsky and Nesvorny (1996) for more details: 
    //     https://ui.adsabs.harvard.edu/abs/1996CeMDA..65..137S/abstract
    int error = reb_frequency_analysis(out, nfreq, -minfreq, minfreq, REB_FREQUENCY_ANALYSIS_FMFT, inp, Nsamples);
    if (error){
        printf("An error occurred during the frequency analysis.\n");
    }
    
    // Output the nfreq most dominate modes
    for (int i=0; i<nfreq; i++){
        double nu = out[0*nfreq+i]*1296000.0/datasep; // frequency in "/year
        double A = out[1*nfreq+i];                    // amplitude error
        double phi = out[2*nfreq+i]/M_PI*180.0;       // phase in deg
        printf("nu = %5.2f\"/yr  A = %0.6f  phi = %5.1f°\n", nu, A, phi);
    }

    // Cleanup
    free(out);
    free(inp);
    reb_simulation_free(r);
}

