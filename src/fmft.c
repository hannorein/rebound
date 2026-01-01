/* This program implements the Frequency Modified Fourier Transform 
   (Sidlichovsky and Nesvorny 1996, Cel. Mech. 65, 137). 
   https://ui.adsabs.harvard.edu/abs/1996CeMDA..65..137S/abstract
   Given a quasi--periodic complex signal X + iY, the algorithm 
   estimates the frequencies (f_j), amplitudes (A_j) and phases 
   (psi_j) in its decomposition:

   X(t) + iY(t) = Sum_j=1^N [ A_j * exp i (f_j * t + psi_j) ] */      

#define FMFT_TOL 1.0e-10 /* MFT NOMINAL PRECISION */
#define FMFT_NEAR 0.     /* MFT OVERLAP EXCLUSION PARAMETER */

#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

#define TWOPI (2.*M_PI)

static void window(double *x, double *y, double *xdata, double *ydata, long ndata);
static void power(double *powsd, double *x, double *y, long ndata);
static void four1(double* data, unsigned long n);
static double bracket(double *powsd, long ndata);
static double golden(double centerf, double width, double *x, double *y, long ndata);
static void phifun(double *xphi, double *yphi, double freq,  double* xdata, double* ydata, long n);
static double phisqr(double freq, double* xdata, double* ydata, long ndata);
static void amph(double *amp, double *phase, double freq, double* xdata, double* ydata, long ndata);
static void sort3(unsigned long n, double* ra, double* rb, double* rc, double* rd);


/* THE MAIN FUNCTION ****************************************************/
void reb_frequency_analysis(double **output, int nfreq, double minfreq, double maxfreq, enum REB_FREQUENCY_ANALYSIS_TYPE type, double *input, long ndata){

    /* 
       In the output array **output: output[i*3], output[i*3+1] 
       and output[i*3+2] are the i-th frequency, amplitude and phase; nfreq is the 
       number of frequencies to be computed (the units are rad/sep, where sep is the 
       `time' separation between i and i+1. The algorithm is  

       Modified Fourier Transform                  if   type = 0;
       Frequency Modified Fourier Transform        if   type = 1;
       FMFT with additional non-linear correction  if   type = 2

       (while the first algorithm is app. 3 times faster than the third one, 
       the third algorithm should be in general much more precise).  
       The computed frequencies are in the range given by minfreq and maxfreq.
       The function returns the number of determined frequencies or 0 in the case
       of error.

       The vector input[j], j = 0 ... 2*ndata-1 (ndata must
       be a power of 2), are the input data X(j) and Y(j).
     */   



    /* ALLOCATION OF VARIABLES */

    double* xdata = malloc(sizeof(double)*ndata);
    double* ydata = malloc(sizeof(double)*ndata);
    double* x = malloc(sizeof(double)*ndata);
    double* y = malloc(sizeof(double)*ndata);
    double* powsd = malloc(sizeof(double)*ndata);

    double* freq = malloc(sizeof(double)*3*(type+1)*nfreq); 
    double* amp = malloc(sizeof(double)*3*(type+1)*nfreq);
    double* phase = malloc(sizeof(double)*3*(type+1)*nfreq);

    double* f = malloc(sizeof(double)*nfreq);
    double* A = malloc(sizeof(double)*nfreq);
    double* psi = malloc(sizeof(double)*nfreq);


    double* Q = malloc(sizeof(double)*nfreq*nfreq);
    double* alpha = malloc(sizeof(double)*nfreq*nfreq);
    double* B = malloc(sizeof(double)*nfreq);


    /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

    for(int l=0; l<=type; l++){
        if(l==0){
            /* SEPARATE REAL AND IMAGINERY PARTS */ 
            for(int j=0;j<ndata;j++){
                xdata[j] = input[j*2];
                ydata[j] = input[j*2+1];
            }
        } else {
            /* GENERATE THE QUASIPERIODIC FUNCTION COMPUTED BY MFT */
            for(int i=0;i<ndata;i++){
                xdata[i] = 0; 
                ydata[i] = 0; 
                for(int k=0;k<nfreq;k++){
                    xdata[i] += amp[(l-1)*nfreq+k]*cos(freq[(l-1)*nfreq+k]*i + phase[(l-1)*nfreq+k]);
                    ydata[i] += amp[(l-1)*nfreq+k]*sin(freq[(l-1)*nfreq+k]*i + phase[(l-1)*nfreq+k]);
                }
            }
        }

        /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
        window(x, y, xdata, ydata, ndata);

        /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
        power(powsd, x, y, ndata);

        double centerf;

        if(l==0){ 
            /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
            while((centerf = bracket(powsd, ndata)) < minfreq || centerf > maxfreq) {
                /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
                f[0] = golden(centerf, TWOPI/ndata, x, y, ndata);

                amph(&A[0], &psi[0], f[0], x, y, ndata);

                for(int j=0;j<ndata;j++){
                    xdata[j] -= A[0]*cos( f[0]*j + psi[0] );
                    ydata[j] -= A[0]*sin( f[0]*j + psi[0] );
                }

                window(x, y, xdata, ydata, ndata);

                power(powsd, x, y, ndata); 
            }   
        }else{ 
            centerf = freq[0];
        }

        /* DETERMINE THE FIRST FREQUENCY */
        f[0] = golden(centerf, TWOPI/ndata, x, y, ndata);

        /* COMPUTE AMPLITUDE AND PHASE */
        amph(&A[0], &psi[0], f[0], x, y, ndata);

        /* SUBSTRACT THE FIRST HARMONIC FROM THE SIGNAL */
        for(int j=0;j<ndata;j++){
            xdata[j] -= A[0]*cos( f[0]*j + psi[0] );
            ydata[j] -= A[0]*sin( f[0]*j + psi[0] );
        }    

        /* HERE STARTS THE MAIN LOOP  *************************************/ 
        Q[0] = 1;
        alpha[0] = 1;

        for(int m=1;m<nfreq;m++){
            /* MULTIPLY SIGNAL BY WINDOW FUNCTION */
            window(x, y, xdata, ydata, ndata);

            /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
            power(powsd, x, y, ndata);

            if(l==0){
                double centerf = bracket(powsd, ndata);
                f[m] = golden(centerf, TWOPI/ndata, x, y, ndata);

                /* CHECK WHETHER THE NEW FREQUENCY IS NOT TOO CLOSE TO ANY PREVIOUSLY
                   DETERMINED ONE */
                int nearfreqflag = 0.;
                for(int k=0;k<m-1;k++){
                    if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata ){
                        nearfreqflag = 1; 
                    }
                }

                /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
                while(f[m] < minfreq || f[m] > maxfreq || nearfreqflag == 1){
                    /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
                    f[m] = golden(centerf, TWOPI/ndata, x, y, ndata);

                    amph(&A[m], &psi[m], f[m], x, y, ndata);

                    for(int j=0;j<ndata;j++){
                        xdata[j] -= A[m]*cos( f[m]*j + psi[m] );
                        ydata[j] -= A[m]*sin( f[m]*j + psi[m] );
                    }

                    /* AND RECOMPUTE THE NEW ONE */
                    window(x, y, xdata, ydata, ndata);

                    power(powsd, x, y, ndata); 

                    centerf = bracket(powsd, ndata); 
                    f[m] = golden(centerf, TWOPI/ndata, x, y, ndata);

                    nearfreqflag = 0.;
                    for(int k=0;k<m-1;k++){
                        if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata ){
                            nearfreqflag = 1; 
                        }
                    }
                }   

            } else {  
                /* DETERMINE THE NEXT FREQUENCY */
                f[m] = golden(freq[m], TWOPI/ndata, x, y, ndata);
            }

            /* COMPUTE ITS AMPLITUDE AND PHASE */
            amph(&A[m], &psi[m], f[m], x, y, ndata);


            /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
            Q[m*nfreq+m] = 1;
            for(int j=0;j<m;j++){
                double fac = (f[m] - f[j]) * (ndata - 1.) / 2.;
                Q[m*nfreq+j] = sin(fac)/fac * M_PI*M_PI / (M_PI*M_PI - fac*fac);
                Q[j*nfreq+m] = Q[m*nfreq+j];
            }

            /* EQUATION (17) */
            for(int k=0;k<m;k++){
                B[k] = 0;
                for(int j=0;j<k;j++)
                    B[k] += -alpha[k*nfreq+j]*Q[m*nfreq+j];
            }

            /* EQUATION (18) */
            alpha[m*nfreq+m] = 1;
            for(int j=0;j<m;j++)
                alpha[m*nfreq+m] -= B[j]*B[j];
            alpha[m*nfreq+m] = 1. / sqrt(alpha[m*nfreq+m]);


            /* EQUATION (19) */
            for(int k=0;k<m;k++){
                alpha[m*nfreq+k] = 0;
                for(int j=k;j<m;j++)
                    alpha[m*nfreq+k] += B[j]*alpha[j*nfreq+k];
                alpha[m*nfreq+k] = alpha[m*nfreq+m]*alpha[m*nfreq+k];
            }

            /* EQUATION (22) */
            for(int i=0;i<ndata;i++){
                double xsum=0; 
                double ysum=0;
                for(int j=0;j<=m;j++){
                    double fac = f[j]*i + (f[m]-f[j])*(ndata-1.)/2. + psi[m];
                    xsum += alpha[m*nfreq+j]*cos(fac);
                    ysum += alpha[m*nfreq+j]*sin(fac);
                }
                xdata[i] -= alpha[m*nfreq+m]*A[m]*xsum;
                ydata[i] -= alpha[m*nfreq+m]*A[m]*ysum;
            }
        }

        /* EQUATION (26) */
        for(int k=0;k<nfreq;k++){
            double xsum=0; 
            double ysum=0;
            for(int j=k;j<nfreq;j++){
                double fac = (f[j]-f[k])*(ndata-1.)/2. + psi[j];
                xsum += alpha[j*nfreq+j]*alpha[j*nfreq+k]*A[j]*cos(fac);
                ysum += alpha[j*nfreq+j]*alpha[j*nfreq+k]*A[j]*sin(fac);
            }
            A[k] = sqrt(xsum*xsum + ysum*ysum);
            psi[k] = atan2(ysum,xsum);
        }

        /* REMEMBER THE COMPUTED VALUES FOR THE FMFT */
        for(int k=0;k<nfreq;k++){
            freq[l*nfreq+k] = f[k];
            amp[l*nfreq+k] = A[k];
            phase[l*nfreq+k] = psi[k];
        }
    }

    /* RETURN THE FINAL FREQUENCIES, AMPLITUDES AND PHASES */ 
    *output = calloc(3*nfreq,sizeof(double));
    switch (type){
        case REB_FREQUENCY_ANALYSIS_MFT:
            for(int k=0;k<nfreq;k++){
                (*output)[0*nfreq+k] = freq[0*nfreq+k];            
                (*output)[1*nfreq+k] = amp[0*nfreq+k];
                (*output)[2*nfreq+k] = phase[0*nfreq+k];
            }
            break;
        case REB_FREQUENCY_ANALYSIS_FMFT:
            for(int k=0;k<nfreq;k++){
                (*output)[0*nfreq+k] = freq[0*nfreq+k] + (freq[0*nfreq+k] - freq[1*nfreq+k]);            
                (*output)[1*nfreq+k] = amp[0*nfreq+k] + (amp[0*nfreq+k] - amp[1*nfreq+k]);
                (*output)[2*nfreq+k] = phase[0*nfreq+k] + (phase[0*nfreq+k] - phase[1*nfreq+k]);
            }
            break;
        case REB_FREQUENCY_ANALYSIS_FMFT2:
            for(int k=0;k<nfreq;k++){
                (*output)[0*nfreq+k] = freq[0*nfreq+k];
                double fac;
                if(fabs((fac = freq[1*nfreq+k] - freq[2*nfreq+k])/freq[1*nfreq+k]) > FMFT_TOL){
                    double tmp = freq[0*nfreq+k] - freq[1*nfreq+k];
                    (*output)[0*nfreq+k] += tmp*tmp / fac;
                }else{ 
                    (*output)[0*nfreq+k] += freq[0*nfreq+k] - freq[1*nfreq+k]; 
                }
                (*output)[1*nfreq+k] = amp[0*nfreq+k];
                if(fabs((fac = amp[1*nfreq+k] - amp[2*nfreq+k])/amp[1*nfreq+k]) > FMFT_TOL){
                    double tmp = amp[0*nfreq+k] - amp[1*nfreq+k];
                    (*output)[1*nfreq+k] += tmp*tmp / fac;
                }else{
                    (*output)[1*nfreq+k] += amp[0*nfreq+k] - amp[1*nfreq+k]; 
                }
                (*output)[2*nfreq+k] = phase[0*nfreq+k];
                if(fabs((fac = phase[1*nfreq+k] - phase[2*nfreq+k])/phase[1*nfreq+k]) > FMFT_TOL){
                    double tmp = phase[0*nfreq+k] - phase[1*nfreq+k];
                    (*output)[2*nfreq+k] += tmp*tmp / fac;
                }else{
                    (*output)[2*nfreq+k] += phase[0*nfreq+k] - phase[1*nfreq+k]; 
                }
            }
            break;
        default:
            printf("REB_FREQUENCY_ANALYSIS_TYPE not implemented.\n");
    }
    for(int k=0;k<nfreq;k++){
        if((*output)[2*nfreq+k] < 0.0) (*output)[2*nfreq+k] += TWOPI;
        if((*output)[2*nfreq+k] >= 2.0*M_PI) (*output)[2*nfreq+k] -= TWOPI;
    }

    // SORT THE FREQUENCIES IN DECREASING ORDER OF AMPLITUDE
    sort3(nfreq, &((*output)[1*nfreq]), &((*output)[0*nfreq]), &((*output)[1*nfreq]), &((*output)[2*nfreq]));

    /* FREE THE ALLOCATED VARIABLES */
    free(xdata);
    free(ydata);
    free(x);
    free(y);
    free(powsd);

    free(freq); 
    free(amp); 
    free(phase); 

    free(f);
    free(A);
    free(psi);

    free(Q); 
    free(alpha);
    free(B);
}

static void window(double *x, double *y, double *xdata, double *ydata, long ndata) {  
    // Hanning window
    for(int j=0;j<ndata;j++) {
        double window = (1. - cos(TWOPI*j / (ndata-1)))*0.5;
        x[j] = xdata[j]*window;
        y[j] = ydata[j]*window;
    }
}


static void power(double *powsd, double *x, double *y, long ndata){
    /* REARRANGES DATA FOR THE FAST FOURIER TRANSFORM, 
       CALLS FFT AND RETURNS POWER SPECTRAL DENSITY */
    double* z = malloc(sizeof(double)*2*ndata);
    for(int j=0;j<ndata;j++){
        z[2*j] = x[j];
        z[2*j+1] = y[j];
    }
    four1(z, ndata);
    for(int j=0;j<ndata;j++){
        powsd[j] = z[2*j]*z[2*j] + z[2*j+1]*z[2*j+1];
    }
    free(z);
}


static void four1(double* data, unsigned long nn){
    /* data[1..2*nn] replaces by DFS, nn must be a power of 2 */
    unsigned long n;

    n=nn<<1;
    unsigned long j=0;
    for(unsigned long i=0;i<n-1;i+=2){ /* bit-reversal section */
        if(j>i){
            double t = data[j];
            data[j] = data[i];
            data[i] = t;
            t = data[j+1];
            data[j+1] = data[i+1];
            data[i+1] = t;
        }
        unsigned long m=n>>1;
        while(m>=2 && j+1>m){
            j-=m;
            m>>=1;
        }
        j+=m;
    }
    /* Danielson-Lanczos section */
    unsigned long mmax=2;
    while(n>mmax){ /* outer ln nn loop */
        unsigned long istep=mmax<<1;
        double theta=TWOPI/mmax; /* initialize */
        double wtemp=sin(0.5*theta);
        double wpr=-2.0*wtemp*wtemp;
        double wpi=sin(theta);
        double wr=1.0;
        double wi=0.0;
        for(unsigned long m=0;m<mmax;m+=2){ /* two inner loops */
            for(int i=m;i<n;i+=istep){
                j=i+mmax; /* D-L formula */
                double tempr=wr*data[j]-wi*data[j+1];
                double tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i]+=tempr;
                data[i+1]+=tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr; /* trig. recurrence */
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

static double bracket(double *powsd, long ndata){
    /* FINDS THE MAXIMUM OF THE POWER SPECTRAL DENSITY  */ 
    // Not sure why this is so complicated. Could just be one loop?
    int maxj = 0;
    double maxpow = 0;

    for(int j=1;j<ndata/2-1;j++){ // Changed end from -2 to -1.
        if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1] && powsd[j] > maxpow){ 
            maxj = j;
            maxpow = powsd[j];
        }
    }

    for(int j=ndata/2+1;j<ndata-1;j++){
        if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1] && powsd[j] > maxpow){ 
            maxj = j;
            maxpow = powsd[j];
        }
    }

    if(powsd[0] > powsd[1] && powsd[0] > powsd[ndata-1] && powsd[0] > maxpow){ 
        maxj = 0;
        maxpow = powsd[0];
    }

    if(maxpow == 0) printf("DFT has no maximum ...");

    if(maxj < ndata/2-1){
        return -TWOPI*maxj / ndata;
    }else{ //  maxj > ndata/2-1
        return -TWOPI*(maxj-ndata) / ndata;
    }
    /* negative signs and TWOPI compensate for the Numerical Recipes 
       definition of the DFT */
}


static double golden(double bx, double width, double* xdata, double* ydata, long n){
    /* calculates the maximum of a function bracketed by ax, bx and cx */
    const double gold_r =  0.6180339887498948482;
    const double gold_c = (1.0 - gold_r);

    double ax = bx-width;
    double cx = bx+width;
    double x0=ax;
    double x3=cx;

    double x1,x2;
    if(fabs(cx-bx) > fabs(bx-ax)){
        x1 = bx;
        x2 = bx + gold_c*(cx-bx);
    } else {
        x2 = bx;
        x1 = bx - gold_c*(bx-ax);
    }

    double f1 = phisqr(x1, xdata, ydata, n);
    double f2 = phisqr(x2, xdata, ydata, n);

    while(fabs(x3-x0) > FMFT_TOL*(fabs(x1)+fabs(x2))){
        if(f2 > f1){
            x0 = x1;
            x1 = x2;
            x2 = gold_r*x1+gold_c*x3;
            f1 = f2;
            f2 = phisqr(x2, xdata, ydata, n);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = gold_r*x2+gold_c*x0;
            f2 = f1;
            f1 = phisqr(x1, xdata, ydata, n);
        }
    }

    if(f1>f2){
        return x1;
    }else{
        return x2;
    }
}

static void amph(double *amp, double *phase, double freq, double* xdata, double* ydata, long ndata){
    /* CALCULATES THE AMPLITUDE AND PHASE */
    double xphi = 0;
    double yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    *amp = sqrt(xphi*xphi + yphi*yphi);
    *phase = atan2(yphi, xphi);
}

static double phisqr(double freq, double* xdata, double* ydata, long ndata){
    /* COMPUTES A SQUARE POWER OF THE FUNCTION PHI */
    double xphi = 0;
    double yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    return xphi*xphi + yphi*yphi;
}

static void phifun(double *xphi, double *yphi, double freq, double* xdata, double* ydata, long n){
    /* COMPUTES THE FUNCTION PHI */   
    double* xdata2 = malloc(sizeof(double)* n);
    double* ydata2 = malloc(sizeof(double)* n);

    xdata2[0] = xdata[0] / 2; ydata2[0] = ydata[0] / 2;
    xdata2[n-1] = xdata[n-1] / 2; ydata2[n-1] = ydata[n-1] / 2;

    for(int i=1;i<n-1;i++){
        xdata2[i] = xdata[i];
        ydata2[i] = ydata[i];
    }

    int nn = n;
    while(nn != 1){
        nn = nn / 2;
        double c = cos(-nn*freq);
        double s = sin(-nn*freq);

        for(int i=0;i<nn;i++){
            int j=i+nn;
            xdata2[i] += c*xdata2[j] - s*ydata2[j];
            ydata2[i] += c*ydata2[j] + s*xdata2[j];
        }
    }

    *xphi = 2*xdata2[0] / (n-1);
    *yphi = 2*ydata2[0] / (n-1);

    free(xdata2);
    free(ydata2);
}

// Workaround because qsort_r is not portable
struct cmp2 {
    unsigned int i;
    double a;
};
static int compare_amp(const void* a, const void* b){
    double aa = ((struct cmp2*)a)->a;
    double ba = ((struct cmp2*)b)->a;
    if (aa>ba) return 1;
    if (aa<ba) return -1;
    return 0;
}

// Sort rb, rc, rd depending on ra
static void sort3(unsigned long n, double* ra, double* rb, double* rc, double* rd){
    struct cmp2* iwksp = malloc(sizeof(struct cmp2)*n);
    double* wksp = malloc(sizeof(double)* n);
    for (unsigned long j=0;j<n;j++){
        iwksp[j].i = j;
        iwksp[j].a = ra[j];
    }
    qsort(iwksp, n, sizeof(struct cmp2), compare_amp);

    for (int j=0;j<n;j++) wksp[j] = rb[j];
    for(int j=0;j<n;j++) rb[j] = wksp[iwksp[n-j-1].i];
    for (int j=0;j<n;j++) wksp[j] = rc[j];
    for(int j=0;j<n;j++) rc[j] = wksp[iwksp[n-j-1].i];
    for (int j=0;j<n;j++) wksp[j] = rd[j];
    for(int j=0;j<n;j++) rd[j] = wksp[iwksp[n-j-1].i];

    free(wksp);
    free(iwksp);
}

