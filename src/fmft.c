/* This program implements the Frequency Modified Fourier Transform 
   (Sidlichovsky and Nesvorny 1997, Cel. Mech. 65, 137). 
   Given a quasi--periodic complex signal X + iY, the algorithm 
   estimates the frequencies (f_j), amplitudes (A_j) and phases 
   (psi_j) in its decomposition:

   X(t) + iY(t) = Sum_j=1^N [ A_j * exp i (f_j * t + psi_j) ] */      

#define FMFT_TOL 1.0e-10 /* MFT NOMINAL PRECISION */
#define FMFT_NEAR 0.     /* MFT OVERLAP EXCLUSION PARAMETER */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TWOPI (2.*M_PI)

static int itemp;
static unsigned long ultemp;
static double dtemp;

#define DSQR(a) ((dtemp=(a)) == 0.0 ? 0.0 : dtemp*dtemp)

#define SHFT3(a,b,c) (a)=(b);(b)=(c)
#define SHFT4(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

#define ULSWAP(a,b) ultemp=(a);(a)=(b);(b)=ultemp

void window(double *x, double *y, double *xdata, double *ydata, long ndata);

void power(double *powsd, double *x, double *y, long ndata);

void four1(double data[], unsigned long n, int isign);

double bracket(double *powsd, long ndata);

double golden(double (*f)(double, double *, double *, long), 
        double leftf, double centerf, double rightf, 
        double *x, double *y, long ndata);

void phifun(double *xphi, double *yphi, double freq,  
        double xdata[], double ydata[], long n);

double phisqr(double freq, double xdata[], double ydata[], long ndata);

void amph(double *amp, double *phase, double freq, 
        double xdata[], double ydata[], long ndata);

void dsort(unsigned long n, double ra[], double rb[], double rc[], double rd[]);

void dindex(unsigned long n, double arr[], unsigned long indx[]);



/* THE MAIN FUNCTION ****************************************************/


int fmft(double **output, int nfreq, double minfreq, double maxfreq, int flag, 
        double **input, long ndata)

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

       The vectors input[1][j] and input[2][j], j = 1 ... ndata (ndata must
       be a power of 2), are the input data X(j-1) and Y(j-1).
     */   

{
    int nearfreqflag;
    long i,j,k,l,m;
    double *powsd;
    double *xdata, *ydata, *x, *y;
    double centerf, leftf, rightf, fac, xsum, ysum;
    double **freq, **amp, **phase, *f, *A, *psi;
    double **Q, **alpha, *B;


    /* ALLOCATION OF VARIABLES */

    xdata = malloc(sizeof(double)*ndata);
    ydata = malloc(sizeof(double)*ndata);
    x = malloc(sizeof(double)*ndata);
    y = malloc(sizeof(double)*ndata);
    powsd = malloc(sizeof(double)* ndata);

    freq = dmatrix(1, 3*flag, 1, nfreq); 
    amp = dmatrix(1, 3*flag, 1, nfreq);
    phase = dmatrix(1, 3*flag, 1, nfreq);

    f = malloc(sizeof(double)* nfreq);
    A = malloc(sizeof(double)* nfreq);
    psi = malloc(sizeof(double)* nfreq);


    Q = dmatrix(1, nfreq, 1, nfreq); 
    alpha = dmatrix(1, nfreq, 1, nfreq);
    B = malloc(sizeof(double)* nfreq);


    /* 1 LOOP FOR MFT, 2 LOOPS FOR FMFT, 3 LOOPS FOR NON-LINEAR FMFT */

    for(l=0; l<flag; l++){

        if(l==0){

            /* SEPARATE REAL AND IMAGINERY PARTS */ 
            for(j=0;j<ndata;j++){
                xdata[j] = input[0][j];
                ydata[j] = input[1][j];
            }

        } else {

            /* GENERATE THE QUASIPERIODIC FUNCTION COMPUTED BY MFT */
            for(i=0;i<ndata;i++){
                xdata[i] = 0; 
                ydata[i] = 0; 
                for(k=0;k<nfreq;k++){
                    xdata[i] += amp[l][k]*cos(freq[l][k]*i + phase[l][k]);
                    ydata[i] += amp[l][k]*sin(freq[l][k]*i + phase[l][k]);
                }
            }

        }

        /* MULTIPLY THE SIGNAL BY A WINDOW FUNCTION, STORE RESULT IN x AND y */
        window(x, y, xdata, ydata, ndata);

        /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
        power(powsd, x, y, ndata);


        if(l==0) 

            /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
            while((centerf = bracket(powsd, ndata)) < minfreq || centerf > maxfreq) {


                /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;

                f[0] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

                amph(&A[0], &psi[0], f[0], x, y, ndata);

                for(j=0;j<ndata;j++){
                    xdata[j] -= A[0]*cos( f[0]*j + psi[0] );
                    ydata[j] -= A[0]*sin( f[0]*j + psi[0] );
                }

                window(x, y, xdata, ydata, ndata);

                power(powsd, x, y, ndata); 
            }   

        else 
            centerf = freq[0][0];

        leftf = centerf - TWOPI / ndata;
        rightf = centerf + TWOPI / ndata;

        /* DETERMINE THE FIRST FREQUENCY */
        f[0] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

        /* COMPUTE AMPLITUDE AND PHASE */
        amph(&A[0], &psi[0], f[0], x, y, ndata);

        /* SUBSTRACT THE FIRST HARMONIC FROM THE SIGNAL */
        for(j=0;j<ndata;j++){
            xdata[j] -= A[0]*cos( f[0]*j + psi[0] );
            ydata[j] -= A[0]*sin( f[0]*j + psi[0] );
        }    

        /* HERE STARTS THE MAIN LOOP  *************************************/ 

        Q[0][0] = 1;
        alpha[0][0] = 1;

        for(m=1;m<nfreq;m++){

            /* MULTIPLY SIGNAL BY WINDOW FUNCTION */
            window(x, y, xdata, ydata, ndata);

            /* COMPUTE POWER SPECTRAL DENSITY USING FAST FOURIER TRANSFORM */
            power(powsd, x, y, ndata);

            if(l==0){

                centerf = bracket(powsd, ndata);

                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;

                f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

                /* CHECK WHETHER THE NEW FREQUENCY IS NOT TOO CLOSE TO ANY PREVIOUSLY
                   DETERMINED ONE */
                nearfreqflag = 0.;
                for(k=0;k<m-1;k++)
                    if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = 1; 

                /* CHECK IF THE FREQUENCY IS IN THE REQUIRED RANGE */
                while(f[m] < minfreq || f[m] > maxfreq || nearfreqflag == 1){

                    /* IF NO, SUBSTRACT IT FROM THE SIGNAL */
                    leftf = centerf - TWOPI / ndata;
                    rightf = centerf + TWOPI / ndata;

                    f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

                    amph(&A[m], &psi[m], f[m], x, y, ndata);

                    for(j=0;j<ndata;j++){
                        xdata[j] -= A[m]*cos( f[m]*j + psi[m] );
                        ydata[j] -= A[m]*sin( f[m]*j + psi[m] );
                    }

                    /* AND RECOMPUTE THE NEW ONE */
                    window(x, y, xdata, ydata, ndata);

                    power(powsd, x, y, ndata); 

                    centerf = bracket(powsd, ndata); 

                    leftf = centerf - TWOPI / ndata;
                    rightf = centerf + TWOPI / ndata;

                    f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

                    nearfreqflag = 0.;
                    for(k=0;k<m-1;k++)
                        if( fabs(f[m] - f[k]) < FMFT_NEAR*TWOPI/ndata )   nearfreqflag = 1; 

                }   

            } else {  

                centerf = freq[0][m];

                leftf = centerf - TWOPI / ndata;
                rightf = centerf + TWOPI / ndata;

                /* DETERMINE THE NEXT FREQUENCY */
                f[m] = golden(phisqr, leftf, centerf, rightf, x, y, ndata);

            }

            /* COMPUTE ITS AMPLITUDE AND PHASE */
            amph(&A[m], &psi[m], f[m], x, y, ndata);


            /* EQUATION (3) in Sidlichovsky and Nesvorny (1997) */
            Q[m][m] = 1;
            for(j=0;j<m-1;j++){
                fac = (f[m] - f[j]) * (ndata - 1.) / 2.;
                Q[m][j] = sin(fac)/fac * M_PI*M_PI / (M_PI*M_PI - fac*fac);
                Q[j][m] = Q[m][j];
            }

            /* EQUATION (17) */
            for(k=0;k<m-1;k++){
                B[k] = 0;
                for(j=0;j<k;j++)
                    B[k] += -alpha[k][j]*Q[m][j];
            }

            /* EQUATION (18) */
            alpha[m][m] = 1;
            for(j=0;j<m-1;j++)
                alpha[m][m] -= B[j]*B[j];
            alpha[m][m] = 1. / sqrt(alpha[m][m]);


            /* EQUATION (19) */
            for(k=0;k<m-1;k++){
                alpha[m][k] = 0;
                for(j=k;j<m-1;j++)
                    alpha[m][k] += B[j]*alpha[j][k];
                alpha[m][k] = alpha[m][m]*alpha[m][k];
            }

            /* EQUATION (22) */
            for(i=0;i<ndata;i++){
                xsum=0; ysum=0;
                for(j=0;j<m;j++){
                    fac = f[j]*i + (f[m]-f[j])*(ndata-1.)/2. + psi[m];
                    xsum += alpha[m][j]*cos(fac);
                    ysum += alpha[m][j]*sin(fac);
                }
                xdata[i] -= alpha[m][m]*A[m]*xsum;
                ydata[i] -= alpha[m][m]*A[m]*ysum;
            }
        }

        /* EQUATION (26) */
        for(k=0;k<nfreq;k++){
            xsum=0; ysum=0;
            for(j=k;j<nfreq;j++){
                fac = (f[j]-f[k])*(ndata-1.)/2. + psi[j];
                xsum += alpha[j][j]*alpha[j][k]*A[j]*cos(fac);
                ysum += alpha[j][j]*alpha[j][k]*A[j]*sin(fac);
            }
            A[k] = sqrt(xsum*xsum + ysum*ysum);
            psi[k] = atan2(ysum,xsum);
        }

        /* REMEMBER THE COMPUTED VALUES FOR THE FMFT */
        for(k=0;k<nfreq;k++){
            freq[l][k] = f[k];
            amp[l][k] = A[k];
            phase[l][k] = psi[k];
        }
    }

    /* RETURN THE FINAL FREQUENCIES, AMPLITUDES AND PHASES */ 

    for(k=0;k<nfreq;k++){
        output[0][k] = freq[0][k];            
        output[1][k] = amp[0][k];
        output[2][k] = phase[0][k];

        if(output[2][k] < -M_PI) output[2][k] += TWOPI;
        if(output[2][k] >= M_PI) output[2][k] -= TWOPI;
    }

    if(flag==1 || flag==2)
        for(k=0;k<nfreq;k++){
            output[3][k] = freq[0][k] + (freq[0][k] - freq[1][k]);            
            output[4][k] = amp[0][k] + (amp[0][k] - amp[1][k]);
            output[5][k] = phase[0][k] + (phase[0][k] - phase[1][k]);

            if(output[5][k] < -M_PI) output[5][k] += TWOPI;
            if(output[5][k] >= M_PI) output[5][k] -= TWOPI;
        }

    if(flag==2)
        for(k=0;k<nfreq;k++){

            output[6][k] = freq[0][k];
            if(fabs((fac = freq[1][k] - freq[2][k])/freq[1][k]) > FMFT_TOL)
                output[6][k] += DSQR(freq[0][k] - freq[1][k]) / fac;
            else 
                output[6][k] += freq[0][k] - freq[1][k]; 

            output[7][k] = amp[0][k];
            if(fabs((fac = amp[1][k] - amp[2][k])/amp[1][k]) > FMFT_TOL)
                output[7][k] += DSQR(amp[0][k] - amp[1][k]) / fac;
            else
                output[7][k] += amp[0][k] - amp[1][k]; 

            output[8][k] = phase[0][k];
            if(fabs((fac = phase[1][k] - phase[2][k])/phase[1][k]) > FMFT_TOL)
                output[8][k] += DSQR(phase[0][k] - phase[1][k]) / fac;
            else
                output[8][k] += phase[0][k] - phase[1][k]; 

            if(output[8][k] < -M_PI) output[8][k] += TWOPI;
            if(output[8][k] >= M_PI) output[8][k] -= TWOPI;
        }

    /* SORT THE FREQUENCIES IN DECREASING ORDER OF AMPLITUDE */
    if(flag==0) 
        dsort(nfreq, output[1], output[0], output[1], output[2]);

    if(flag==1){
        dsort(nfreq, output[4], output[0], output[1], output[2]);
        dsort(nfreq, output[4], output[3], output[4], output[5]);
    }

    if(flag==2){
        dsort(nfreq, output[7], output[0], output[1], output[2]);
        dsort(nfreq, output[7], output[3], output[4], output[5]);   
        dsort(nfreq, output[7], output[6], output[7], output[8]);
    }

    /* FREE THE ALLOCATED VARIABLES */
    free(xdata);
    free(ydata);
    free(x);
    free(y);
    free(powsd);

    free_dmatrix(freq, 1, 3*flag, 1, nfreq); 
    free_dmatrix(amp, 1, 3*flag, 1, nfreq);
    free_dmatrix(phase, 1, 3*flag, 1, nfreq);

    free(f);
    free(A);
    free(psi);

    free_dmatrix(Q, 1, nfreq, 1, nfreq); 
    free_dmatrix(alpha, 1, nfreq, 1, nfreq);
    free(B);

    return 1;
}


void window(double *x, double *y, double *xdata, double *ydata, long ndata)

    /* MULTIPLIES DATA BY A WINDOW FUNCTION */      
{  
    long j;
    double window;

    for(j=0;j<ndata;j++) {

        window = TWOPI*j / (ndata-1);
        window = (1. - cos(window)) / 2.;

        x[j] = xdata[j]*window;
        y[j] = ydata[j]*window;

    }
}


void power(double *powsd, double *x, double *y, long ndata)

    /* REARRANGES DATA FOR THE FAST FOURIER TRANSFORM, 
       CALLS FFT AND RETURNS POWER SPECTRAL DENSITY */

{
    long j;
    double *z;

    z = malloc(sizeof(double)*2*ndata);

    for(j=0;j<ndata;j++){
        z[2*j] = x[j];
        z[2*j+1] = y[j];
    }

    four1(z, ndata, 1);

    for(j=0;j<ndata;j++)
        powsd[j] = DSQR(z[2*j]) + DSQR(z[2*j+1]);

    free(z);
}


void four1(double data[], unsigned long nn, int isign)

    /* data[1..2*nn] replaces by DFS, nn must be a power of 2 */

{
    unsigned long n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta; /* double for recurrences */
    double tempr,tempi;

    n=nn<<1;
    j=0;
    for(i=0;i<n-1;i+=2){ /* bit-reversal section */
        if(j>i){
            double t = data[j];
            data[j] = data[i];
            data[i] = t;
            t = data[j+1];
            data[j+1] = data[i+1];
            data[i+1] = t;
        }
        m=n>>1;
        while(m>=2 && j>m){
            j-=m;
            m>>=1;
        }
        j+=m;
    }
    /* Danielson-Lanczos section */
    mmax=2;
    while(n>mmax){ /* outer ln nn loop */
        istep=mmax<<1;
        theta=isign*(TWOPI/mmax); /* initialize */
        wtemp=sin(0.5*theta);
        wpr=-2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for(m=0;m<mmax-1;m+=2){ /* two inner loops */
            for(i=m;i<n;i+=istep){
                j=i+mmax; /* D-L formula */
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
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

double bracket(double *powsd, long ndata)

    /* FINDS THE MAXIMUM OF THE POWER SPECTRAL DENSITY  */ 

{
    long j, maxj;
    double freq, maxpow;

    maxj = 0;
    maxpow = 0;

    for(j=1;j<ndata/2-2;j++)
        if(powsd[j] > powsd[j-1] && powsd[j] > powsd[j+1])
            if(powsd[j] > maxpow){ 
                maxj = j;
                maxpow = powsd[j];
            }  

    for(j=ndata/2+1;j<ndata-1;j++)
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

    if(maxj < ndata/2) freq = -(maxj-1);  
    if(maxj > ndata/2) freq = -(maxj-ndata-1);

    return (TWOPI*freq / ndata);

    /* negative signs and TWOPI compensate for the Numerical Recipes 
       definition of the DFT */
}

#define GOLD_R 0.61803399
#define GOLD_C (1.0 - GOLD_R)

double golden(double (*f)(double, double *, double *, long), 
        double ax, double bx, double cx,
        double xdata[], double ydata[], long n)

    /* calculates the maximum of a function bracketed by ax, bx and cx */

{
    double f1,f2,x0,x1,x2,x3;

    x0=ax;
    x3=cx;

    if(fabs(cx-bx) > fabs(bx-ax)){
        x1 = bx;
        x2 = bx + GOLD_C*(cx-bx);
    } else {
        x2 = bx;
        x1 = bx - GOLD_C*(bx-ax);
    }

    f1 = (*f)(x1, xdata, ydata, n);
    f2 = (*f)(x2, xdata, ydata, n);

    while(fabs(x3-x0) > FMFT_TOL*(fabs(x1)+fabs(x2))){
        if(f2 > f1){
            SHFT4(x0,x1,x2,GOLD_R*x1+GOLD_C*x3);
            SHFT3(f1,f2,(*f)(x2, xdata, ydata, n));
        } else {
            SHFT4(x3,x2,x1,GOLD_R*x2+GOLD_C*x0);
            SHFT3(f2,f1,(*f)(x1, xdata, ydata, n));
        }
    }

    if(f1>f2) return x1;
    else return x2;
}

void amph(double *amp, double *phase, double freq, 
        double xdata[], double ydata[], long ndata){

    /* CALCULATES THE AMPLITUDE AND PHASE */

    double xphi,yphi;

    xphi = yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    *amp = sqrt(xphi*xphi + yphi*yphi);
    *phase = atan2(yphi, xphi);
}

double phisqr(double freq, double xdata[], double ydata[], long ndata)

    /* COMPUTES A SQUARE POWER OF THE FUNCTION PHI */

{	
    double xphi,yphi;

    xphi = yphi = 0;

    phifun(&xphi, &yphi, freq, xdata, ydata, ndata);

    return xphi*xphi + yphi*yphi;
}

void phifun(double *xphi, double *yphi, double freq,  
        double xdata[], double ydata[], long n)

    /* COMPUTES THE FUNCTION PHI */   

{
    long i, j, nn;
    double c, s, *xdata2, *ydata2;

    xdata2 = malloc(sizeof(double)* n);
    ydata2 = malloc(sizeof(double)* n);

    xdata2[0] = xdata[0] / 2; ydata2[0] = ydata[0] / 2;
    xdata2[n-1] = xdata[n-1] / 2; ydata2[n-1] = ydata[n-1] / 2;

    for(i=0;i<n-2;i++){
        xdata2[i] = xdata[i];
        ydata2[i] = ydata[i];
    }

    nn = n;

    while(nn != 1){

        nn = nn / 2;

        c = cos(-nn*freq);
        s = sin(-nn*freq);

        for(i=0;i<nn;i++){
            j=i+nn;
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

void dsort(unsigned long n, double ra[], double rb[], double rc[], double rd[])

    /* SORTING PROCEDURE FROM NUMERICAL RECIPES */

{
    unsigned long j,*iwksp,n2;
    double *wksp;

    n2 = n+1;
    iwksp = malloc(sizeof(unsigned long)*n);
    wksp = malloc(sizeof(double)* n);

    dindex(n, ra, iwksp);

    for (j=0;j<n;j++) wksp[j] = rb[j];
    for(j=0;j<n;j++) rb[j] = wksp[iwksp[n2-j]];
    for (j=0;j<n;j++) wksp[j] = rc[j];
    for(j=0;j<n;j++) rc[j] = wksp[iwksp[n2-j]];
    for (j=0;j<n;j++) wksp[j] = rd[j];
    for(j=0;j<n;j++) rd[j] = wksp[iwksp[n2-j]];

    free(wksp);
    free(iwksp);
}


void dindex(unsigned long n, double arr[], unsigned long indx[])
{
    unsigned long i,indxt,ir=n,j,k,l=0;
    int jstack=0,*istack;
    double a;

    istack=malloc(sizeof(int)*SORT_NSTACK);
    for (j=0;j<n;j++) indx[j]=j;
    for(;;){
        if(ir-l < SORT_M) {
            for(j=l+1;j<ir;j++) {
                indxt=indx[j];
                a=arr[indxt];
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
            k=(l+ir) >> 1;
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
            a=arr[indxt];
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






