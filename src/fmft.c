/* This program implements the Frequency Modified Fourier Transform 
   (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137). 
   Given a quasi--periodic complex signal X + iY, the algorithm 
   estimates the frequencies (f_j), amplitudes (A_j) and phases 
   (psi_j) in its decomposition:

   X(t) + iY(t) = Sum_j=1^N [ A_j * exp i (f_j * t + psi_j) ] */      

#define FMFT_TOL 1.0e-10 /* MFT NOMINAL PRECISION */
#define FMFT_NEAR 0.     /* MFT OVERLAP EXCLUSION PARAMETER */

#include "rebound.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TWOPI (2.*M_PI)
#define SHFT3(a,b,c) (a)=(b);(b)=(c)
#define SHFT4(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
#define ULSWAP(a,b) ultemp=(a);(a)=(b);(b)=ultemp

static void window(double *x, double *y, double *xdata, double *ydata, long ndata);
static void power(double *powsd, double *x, double *y, long ndata);
static void four1(double data[], unsigned long n);
static double bracket(double *powsd, long ndata);
static double golden(double leftf, double centerf, double rightf, double *x, double *y, long ndata);
static void phifun(double *xphi, double *yphi, double freq,  double xdata[], double ydata[], long n);
static double phisqr(double freq, double xdata[], double ydata[], long ndata);
static void amph(double *amp, double *phase, double freq, double xdata[], double ydata[], long ndata);
static void dsort(unsigned long n, double ra[], double rb[], double rc[], double rd[]);
static void dindex(unsigned long n, double arr[], unsigned long indx[]);


/* THE MAIN FUNCTION ****************************************************/
int reb_fmft(double *output, int nfreq, double minfreq, double maxfreq, int flag, double *input, long ndata){

    /* 
       In the output array **output: output[3*flag-2][i], output[3*flag-1][i] 
       and output[3*flag][i] are the i-th frequency, amplitude and phase; nfreq is the 
       number of frequencies to be computed (the units are rad/sep, where sep is the 
       `time' separation between i and i+1. The algorithm is  

       Modified Fourier Transform                  if   flag = 0;
       Frequency Modified Fourier Transform        if   flag = 1;
       FMFT with additional non-linear correction  if   flag = 2

       (while the first algorithm is app. 3 times faster than the third one, 
       the third algorithm should be in general much more precise).  
       The computed frequencies are in the range given by minfreq and maxfreq.
       The function returns the number of determined frequencies or 0 in the case
       of error.

       The vector input[j], j = 0 ... 2*ndata-1 (ndata must
       be a power of 2), are the input data X(j) and Y(j).
     */   

    double centerf, leftf, rightf;


    /* ALLOCATION OF VARIABLES */

    double* xdata = malloc(sizeof(double)*ndata);
    double* ydata = malloc(sizeof(double)*ndata);
    double* x = malloc(sizeof(double)*ndata);
    double* y = malloc(sizeof(double)*ndata);
    double* powsd = malloc(sizeof(double)* ndata);

    double* freq = malloc(sizeof(double)*3*(flag+1)*nfreq); 
    double* amp = malloc(sizeof(double)*3*(flag+1)*nfreq);
    double* phase = malloc(sizeof(double)*3*(flag+1)*nfreq);

    double* f = malloc(sizeof(double)* nfreq);
    double* A = malloc(sizeof(double)* nfreq);
    double* psi = malloc(sizeof(double)* nfreq);


    double* Q = malloc(sizeof(double)*nfreq*nfreq);
    double* alpha = malloc(sizeof(double)*nfreq*nfreq);
    double* B = malloc(sizeof(double)* nfreq);


    /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

    for(int l=0; l<=flag; l++){
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
                    xdata[i] += amp[l*nfreq+k]*cos(freq[l*nfreq+k]*i + phase[l*nfreq+k]);
                    ydata[i] += amp[l*nfreq+k]*sin(freq[l*nfreq+k]*i + phase[l*nfreq+k]);
                }
            }
        }

        /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
        window(x, y, xdata, ydata, ndata);

        /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
        power(powsd, x, y, ndata);

        if(l==0){ 
            /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
            while((centerf = bracket(powsd, ndata)) < minfreq || centerf > maxfreq) {
                /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;

                f[0] = golden(leftf, centerf, rightf, x, y, ndata);

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
        leftf = centerf - TWOPI / ndata;
        rightf = centerf + TWOPI / ndata;

        /* DETERMINE THE FIRST FREQUENCY */
        f[0] = golden(leftf, centerf, rightf, x, y, ndata);

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
                centerf = bracket(powsd, ndata);
                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;

                f[m] = golden(leftf, centerf, rightf, x, y, ndata);

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
                    leftf = centerf - TWOPI / ndata;
                    rightf = centerf + TWOPI / ndata;

                    f[m] = golden(leftf, centerf, rightf, x, y, ndata);

                    amph(&A[m], &psi[m], f[m], x, y, ndata);

                    for(int j=0;j<ndata;j++){
                        xdata[j] -= A[m]*cos( f[m]*j + psi[m] );
                        ydata[j] -= A[m]*sin( f[m]*j + psi[m] );
                    }

                    /* AND RECOMPUTE THE NEW ONE */
                    window(x, y, xdata, ydata, ndata);

                    power(powsd, x, y, ndata); 

                    centerf = bracket(powsd, ndata); 

                    leftf = centerf - TWOPI / ndata;
                    rightf = centerf + TWOPI / ndata;

                    f[m] = golden(leftf, centerf, rightf, x, y, ndata);

                    nearfreqflag = 0.;
                    for(int k=0;k<m-1;k++){
                        if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata ){
                            nearfreqflag = 1; 
                        }
                    }
                }   

            } else {  
                centerf = freq[m];
                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;
                /* DETERMINE THE NEXT FREQUENCY */
                f[m] = golden(leftf, centerf, rightf, x, y, ndata);
            }

            /* COMPUTE ITS AMPLITUDE AND PHASE */
            amph(&A[m], &psi[m], f[m], x, y, ndata);


            /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
            Q[m*nfreq+m] = 1;
            for(int j=0;j<m-1;j++){
                double fac = (f[m] - f[j]) * (ndata - 1.) / 2.;
                Q[m*nfreq+j] = sin(fac)/fac * M_PI*M_PI / (M_PI*M_PI - fac*fac);
                Q[j*nfreq+m] = Q[m*nfreq+j];
            }

            /* EQUATION (17) */
            for(int k=0;k<m-1;k++){
                B[k] = 0;
                for(int j=0;j<k;j++)
                    B[k] += -alpha[k*nfreq+j]*Q[m*nfreq+j];
            }

            /* EQUATION (18) */
            alpha[m*nfreq+m] = 1;
            for(int j=0;j<m-1;j++)
                alpha[m*nfreq+m] -= B[j]*B[j];
            alpha[m*nfreq+m] = 1. / sqrt(alpha[m*nfreq+m]);


            /* EQUATION (19) */
            for(int k=0;k<m-1;k++){
                alpha[m*nfreq+k] = 0;
                for(int j=k;j<m-1;j++)
                    alpha[m*nfreq+k] += B[j]*alpha[j*nfreq+k];
                alpha[m*nfreq+k] = alpha[m*nfreq+m]*alpha[m*nfreq+k];
            }

            /* EQUATION (22) */
            for(int i=0;i<ndata;i++){
                double xsum=0; 
                double ysum=0;
                for(int j=0;j<m;j++){
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

    for(int k=0;k<nfreq;k++){
        output[0*nfreq+k] = freq[k];            
        output[1*nfreq+k] = amp[0*nfreq+k];
        output[2*nfreq+k] = phase[0*nfreq+k];

        if(output[2*nfreq+k] < -M_PI) output[2*nfreq+k] += TWOPI;
        if(output[2*nfreq+k] >= M_PI) output[2*nfreq+k] -= TWOPI;
    }

    if(flag==1 || flag==2){
        for(int k=0;k<nfreq;k++){
            output[3*nfreq+k] = freq[k] + (freq[k] - freq[1*nfreq+k]);            
            output[4*nfreq+k] = amp[0*nfreq+k] + (amp[0*nfreq+k] - amp[1*nfreq+k]);
            output[5*nfreq+k] = phase[0*nfreq+k] + (phase[0*nfreq+k] - phase[1*nfreq+k]);

            if(output[5*nfreq+k] < -M_PI) output[5*nfreq+k] += TWOPI;
            if(output[5*nfreq+k] >= M_PI) output[5*nfreq+k] -= TWOPI;
        }
    }
    if(flag==2){
        for(int k=0;k<nfreq;k++){

            output[6*nfreq+k] = freq[0*nfreq+k];
            double fac;
            if(fabs((fac = freq[1*nfreq+k] - freq[2*nfreq+k])/freq[1*nfreq+k]) > FMFT_TOL){
                double tmp = freq[0*nfreq+k] - freq[1*nfreq+k];
                output[6*nfreq+k] += tmp*tmp / fac;
            }else 
                output[6*nfreq+k] += freq[0*nfreq+k] - freq[1*nfreq+k]; 

            output[7*nfreq+k] = amp[0*nfreq+k];
            if(fabs((fac = amp[1*nfreq+k] - amp[2*nfreq+k])/amp[1*nfreq+k]) > FMFT_TOL){
                double tmp = amp[0*nfreq+k] - amp[1*nfreq+k];
                output[7*nfreq+k] += tmp*tmp / fac;
            }else
                output[7*nfreq+k] += amp[0*nfreq+k] - amp[1*nfreq+k]; 

            output[8*nfreq+k] = phase[0*nfreq+k];
            if(fabs((fac = phase[1*nfreq+k] - phase[2*nfreq+k])/phase[1*nfreq+k]) > FMFT_TOL){
                double tmp = phase[0*nfreq+k] - phase[1*nfreq+k];
                output[8*nfreq+k] += tmp*tmp / fac;
            }else
                output[8*nfreq+k] += phase[0*nfreq+k] - phase[1*nfreq+k]; 

            if(output[8*nfreq+k] < -M_PI) output[8*nfreq+k] += TWOPI;
            if(output[8*nfreq+k] >= M_PI) output[8*nfreq+k] -= TWOPI;
        }
    }

   // /* SORT THE FREQUENCIES IN DECREASING ORDER OF AMPLITUDE */
    if(flag==0) 
        dsort(nfreq, &output[1*nfreq], &output[0*nfreq], &output[1*nfreq], &output[2*nfreq]);

    if(flag==1){
        dsort(nfreq, &output[4*nfreq], &output[0*nfreq], &output[1*nfreq], &output[2*nfreq]);
        dsort(nfreq, &output[4*nfreq], &output[3*nfreq], &output[4*nfreq], &output[5*nfreq]);
    }

    if(flag==2){
        dsort(nfreq, &output[7*nfreq], &output[0*nfreq], &output[1*nfreq], &output[2*nfreq]);
        dsort(nfreq, &output[7*nfreq], &output[3*nfreq], &output[4*nfreq], &output[5*nfreq]);   
        dsort(nfreq, &output[7*nfreq], &output[6*nfreq], &output[7*nfreq], &output[8*nfreq]);
    }

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

    return 1;
}


void window(double *x, double *y, double *xdata, double *ydata, long ndata) {  
    /* MULTIPLIES DATA BY A WINDOW FUNCTION */      
    for(int j=0;j<ndata;j++) {
        double window = TWOPI*j / (ndata-1);
        window = (1. - cos(window)) / 2.;
        x[j] = xdata[j]*window;
        y[j] = ydata[j]*window;
    }
}


void power(double *powsd, double *x, double *y, long ndata){
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


void four1(double* data, unsigned long nn){
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

double bracket(double *powsd, long ndata){
    /* FINDS THE MAXIMUM OF THE POWER SPECTRAL DENSITY  */ 

    int maxj = 0;
    double maxpow = 0;

    for(int j=1;j<ndata/2-2;j++)
        if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1])
            if(powsd[j] > maxpow){ 
                maxj = j;
                maxpow = powsd[j];
            }  

    for(int j=ndata/2+1;j<ndata-1;j++)
        if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1])
            if(powsd[j] > maxpow){ 
                maxj = j;
                maxpow = powsd[j];
            }  

    if(powsd[1] > powsd[2] && powsd[1] > powsd[ndata])
        if(powsd[1] > maxpow){ 
            maxj = 1;
            maxpow = powsd[1];
        }

    if(maxpow == 0) printf("DFT has no maximum ...");

    double freq = -1.11111; // fix this
    if(maxj < ndata/2-1) freq = -(maxj);  
    if(maxj > ndata/2-1) freq = -(maxj-ndata);
    if (freq==-1.11111) printf("Error occured -1.11111\n");

    return (TWOPI*freq / ndata);

    /* negative signs and TWOPI compensate for the Numerical Recipes 
       definition of the DFT */
}

#define GOLD_R 0.61803399
#define GOLD_C (1.0 - GOLD_R)

double golden(double ax, double bx, double cx, double xdata[], double ydata[], long n){
    /* calculates the maximum of a function bracketed by ax, bx and cx */

    double x0=ax;
    double x3=cx;

    double x1,x2;
    if(fabs(cx-bx) > fabs(bx-ax)){
        x1 = bx;
        x2 = bx + GOLD_C*(cx-bx);
    } else {
        x2 = bx;
        x1 = bx - GOLD_C*(bx-ax);
    }

    double f1 = phisqr(x1, xdata, ydata, n);
    double f2 = phisqr(x2, xdata, ydata, n);

    while(fabs(x3-x0) > FMFT_TOL*(fabs(x1)+fabs(x2))){
        if(f2 > f1){
            SHFT4(x0,x1,x2,GOLD_R*x1+GOLD_C*x3);
            SHFT3(f1,f2,phisqr(x2, xdata, ydata, n));
        } else {
            SHFT4(x3,x2,x1,GOLD_R*x2+GOLD_C*x0);
            SHFT3(f2,f1,phisqr(x1, xdata, ydata, n));
        }
    }

    if(f1>f2) return x1;
    else return x2;
}

void amph(double *amp, double *phase, double freq, double xdata[], double ydata[], long ndata){
    /* CALCULATES THE AMPLITUDE AND PHASE */
    double xphi = 0;
    double yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    *amp = sqrt(xphi*xphi + yphi*yphi);
    *phase = atan2(yphi, xphi);
}

double phisqr(double freq, double xdata[], double ydata[], long ndata){
    /* COMPUTES A SQUARE POWER OF THE FUNCTION PHI */
    double xphi = 0;
    double yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    return xphi*xphi + yphi*yphi;
}

void phifun(double *xphi, double *yphi, double freq, double xdata[], double ydata[], long n){
    /* COMPUTES THE FUNCTION PHI */   
    double* xdata2 = malloc(sizeof(double)* n);
    double* ydata2 = malloc(sizeof(double)* n);

    xdata2[0] = xdata[0] / 2; ydata2[0] = ydata[0] / 2;
    xdata2[n-1] = xdata[n-1] / 2; ydata2[n-1] = ydata[n-1] / 2;

    for(int i=0;i<n-2;i++){
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

#define SORT_M 7 
#define SORT_NSTACK 50

void dsort(unsigned long n, double* ra, double* rb, double* rc, double* rd){
    /* SORTING PROCEDURE FROM NUMERICAL RECIPES */
    unsigned long* iwksp = malloc(sizeof(unsigned long)*n);
    double* wksp = malloc(sizeof(double)* n);

    dindex(n, ra, iwksp);

    for (int j=0;j<n;j++) wksp[j] = rb[j];
    for(int j=0;j<n;j++) rb[j] = wksp[iwksp[n-j-1]];
    for (int j=0;j<n;j++) wksp[j] = rc[j];
    for(int j=0;j<n;j++) rc[j] = wksp[iwksp[n-j-1]];
    for (int j=0;j<n;j++) wksp[j] = rd[j];
    for(int j=0;j<n;j++) rd[j] = wksp[iwksp[n-j-1]];

    free(wksp);
    free(iwksp);
}


void dindex(unsigned long n, double arr[], unsigned long indx[]) {
    unsigned long i,indxt,ir=n,j,l=0;
    int jstack=0;

    int* istack=malloc(sizeof(int)*SORT_NSTACK);
    for (j=0;j<n;j++) indx[j]=j;
    for(;;){
        if(ir-l < SORT_M) {
            for(j=l+1;j<ir;j++) {
                indxt=indx[j];
                double a=arr[indxt];
                for(i=j-1;i>=0;i--) {
                    if(arr[indx[i]] <= a) break;
                    indx[i+1]=indx[i];
                }
                indx[i+1]=indxt;
            }
            if (jstack == 0) break;
            ir=istack[jstack--];
            l=istack[jstack--];
        } else {
            unsigned long k=(l+ir) >> 1;
            unsigned long ultemp;
            ULSWAP(indx[k],indx[l+1]);
            if (arr[indx[l+1]] > arr[indx[ir]]) {
                ULSWAP(indx[l+1],indx[ir]);
            }
            if (arr[indx[l]] > arr[indx[ir]]) {
                ULSWAP(indx[l],indx[ir]);
            }
            if (arr[indx[l+1]] > arr[indx[l]]) {
                ULSWAP(indx[l+1],indx[l]);
            }
            i=l+1;
            j=ir;
            indxt=indx[l];
            double a=arr[indxt];
            for(;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if(j < i) break;
                ULSWAP(indx[i],indx[j]);
            }
            indx[l]=indx[j];
            indx[j]=indxt;
            jstack += 2;
            if (jstack > SORT_NSTACK) printf("SORT_NSTACK too small.");
            if(ir-i+1 >= j-l) {
                istack[jstack]=ir;
                istack[jstack-1]=i;
                ir=j-1;
            } else {
                istack[jstack]=j-1;
                istack[jstack-1]=l;
                l=i;
            }
        }
    }
    free(istack);
}






