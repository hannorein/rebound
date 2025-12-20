#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]) {
    int Nsamples = 1024;
    int nfreq = 4;
    double* input = malloc(sizeof(double)*Nsamples*2);
    double* output = malloc(sizeof(double)*10*nfreq);
    for (int i=0; i<Nsamples; i++){
        input[i*2+0] = sin(0.1*2.0*M_PI*i);
        input[i*2+1] = cos(0.1*2.0*M_PI*i);
    }
    reb_fmft(output, nfreq,0.01,10.,0,input,Nsamples);
    for (int i=0; i<nfreq; i++){
        printf("%d %f %f %f\n", i, output[0*nfreq+i], output[1*nfreq+i], output[2*nfreq+i]);
    }
    free(input);
    free(output);
}

