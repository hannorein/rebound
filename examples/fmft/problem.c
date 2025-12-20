#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]) {
    int Nsamples = 64;
    int nfreq = 4;
    double* input = malloc(sizeof(double)*Nsamples*2);
    double* output = malloc(sizeof(double)*10*nfreq);
    double omega = 10;
    double omega2 = 12;
    double omega3 = 18;
    double datasep = 0.0123;
    for (int i=0; i<Nsamples; i++){
        input[i*2+0] = sin(omega*2.0*M_PI*i*datasep);
        input[i*2+1] = cos(omega*2.0*M_PI*i*datasep);
        input[i*2+0] += 0.1*sin(omega2*2.0*M_PI*i*datasep);
        input[i*2+1] += 0.1*cos(omega2*2.0*M_PI*i*datasep);
        input[i*2+0] += 0.05*sin(omega3*2.0*M_PI*i*datasep);
        input[i*2+1] += 0.05*cos(omega3*2.0*M_PI*i*datasep);
    }
    reb_fmft(output, nfreq,-40,40.,0,input,Nsamples);
    for (int i=0; i<nfreq; i++){
        printf("%d %f %f %f %f %f %f\n", i, output[0*nfreq+i]/datasep/2./M_PI, output[1*nfreq+i], output[2*nfreq+i], output[3*nfreq+i]/datasep/2./M_PI, output[4*nfreq+i], output[5*nfreq+i]);
    }
    free(input);
    free(output);
}

