#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Secular modes for Jupiter. Taken from Laskar (1990).
double nu5[] = {4.2488163, 28.2206942, 3.0895148, 52.1925732, 27.0613982, 29.3799573, 28.8679427, 27.5734578, 5.4070444, 0.6671228}; // frequency, "/yr
double A5[] = {44119.0e-6, 15750.0e-6, 1800.0e-6, 516.0e-6, 183.0e-6, 178.0e-6, 107.0e-6, 95.0e-6, 62.0e-6, 58.0e-6}; // amplitude
double phi5[] = {30.676, 308.112, 121.362, 45.551, 218.696, 217.460, 32.614, 43.733, 116.984, 74.116}; // phase, deg

int main(int argc, char* argv[]) {
    // Create artificial test signal based on Laskar (1990) model.
    int Nsamples = 32768;
    int nfreq = 10;
    double* input = calloc(Nsamples*2, sizeof(double));
    double datasep = 120000.0/365.25*2.0*M_PI; // 120000 days in units of year/2pi
    for (int i=0; i<Nsamples; i++){
        for (int j=0; j<nfreq; j++){
            double nu = nu5[j]/1296000.0; // frequency in units of radians/(year/2pi)
            input[i*2+0] += A5[j]*cos(nu*i*datasep+phi5[j]/180.0*M_PI);
            input[i*2+1] += A5[j]*sin(nu*i*datasep+phi5[j]/180.0*M_PI);
        }
    }
    double* output = malloc(sizeof(double)*10*nfreq);
    reb_fmft(output, nfreq,-40e-2,40.e-2,0,input,Nsamples);
   
    // Check accuracy 
    for (int i=0; i<nfreq; i++){
        assert(fabs(output[0*nfreq+i]/M_PI*648000.0/datasep*2*M_PI-nu5[i]) < 1e-3); // Frequency accuracy better than 1e-3"/year
        assert(fabs((output[1*nfreq+i]-A5[i])/A5[i]) < 1e-2); // Amplitude accuracy better than 1%
        double phi = output[2*nfreq+i]/M_PI*180.0;
        if (phi<0) phi+= 360.0;
        assert(fabs(phi-phi5[i])<1.0); // Phase accuracy better than 1 deg
    
    }
    free(input);
    free(output);
}

