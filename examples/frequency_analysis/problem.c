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
    // Check all 3 types of frequency analysis implemented
    for (enum REB_FREQUENCY_ANALYSIS_TYPE type=0;type<3;type++){
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
        double* output;
        double minfreq = 60.0/1296000.0*datasep;
        reb_frequency_analysis(&output, nfreq,-minfreq,minfreq,type,input,Nsamples);

        // Check accuracy 
        double max_nu_error = 0.0;
        double max_A_error = 0.0;
        double max_phi_error = 0.0;
        for (int i=0; i<nfreq; i++){
            double nu_error = fabs(output[0*nfreq+i]*1296000.0/datasep-nu5[i]); // frequency error in "/year
            if (nu_error > max_nu_error) max_nu_error = nu_error;
            double A_error = fabs((output[1*nfreq+i]-A5[i])/A5[i]); // relative amplitude error
            if (A_error > max_A_error) max_A_error = A_error;
            double phi_error = output[2*nfreq+i]/M_PI*180.0 - phi5[i];
            if (phi_error<-180.0) phi_error+= 360.0;
            if (phi_error>180.0) phi_error-= 360.0;
            phi_error = fabs(phi_error);
            if (phi_error > max_phi_error) max_phi_error = phi_error;
        }
        printf("Flag %d\n", type);
        printf("Max frequency error:          %e \"/year\n", max_nu_error);
        printf("Max relative amplitude error: %e\n", max_A_error);
        printf("Max phase error:              %e deg\n", max_phi_error);
        switch (type){
            case REB_FREQUENCY_ANALYSIS_MFT:
                assert(max_nu_error<3e-4);
                assert(max_A_error<2e-3);
                assert(max_phi_error<5e-1);
                break;
            case REB_FREQUENCY_ANALYSIS_FMFT:
                assert(max_nu_error<4e-6);
                assert(max_A_error<1e-5);
                assert(max_phi_error<6e-3);
                break;
            case REB_FREQUENCY_ANALYSIS_FMFT2:
                assert(max_nu_error<2e-8);
                assert(max_A_error<3e-7);
                assert(max_phi_error<3e-5);
                break;
        }
        free(input);
        free(output);
    }
}

