/**
 * @file    integrator_whfast.c
 * @brief   WHFAST integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the WHFast integration scheme.  
 * Described in Rein & Tamayo 2015.
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Daniel Tamayo
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include "rebound.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator.h"
#include "integrator_whfast.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

// Corrector coefficients
const static double reb_whfast_corrector_a_1 = 0.41833001326703777398908601289259374469640768464934;
const static double reb_whfast_corrector_a_2 = 0.83666002653407554797817202578518748939281536929867;
const static double reb_whfast_corrector_a_3 = 1.2549900398011133219672580386777812340892230539480;
const static double reb_whfast_corrector_a_4 = 1.6733200530681510959563440515703749787856307385973;
const static double reb_whfast_corrector_a_5 = 2.0916500663351888699454300644629687234820384232467;
const static double reb_whfast_corrector_b_31 = -0.024900596027799867499350357910273437184309981229127;
const static double reb_whfast_corrector_b_51 = -0.0083001986759332891664501193034244790614366604097090;
const static double reb_whfast_corrector_b_52 = 0.041500993379666445832250596517122395307183302048545;
const static double reb_whfast_corrector_b_71 = 0.0024926811426922105779030593952776964450539008582219;
const static double reb_whfast_corrector_b_72 = -0.018270923246702131478062356884535264841652263842597;
const static double reb_whfast_corrector_b_73 = 0.053964399093127498721765893493510877532452806339655;
const static double reb_whfast_corrector_b_111 = 0.00020361579647854651301632818774633716473696537436847;
const static double reb_whfast_corrector_b_112 = -0.0023487215292295354188307328851055489876255097419754;
const static double reb_whfast_corrector_b_113 = 0.012309078592019946317544564763237909911330686448336;
const static double reb_whfast_corrector_b_114 = -0.038121613681288650508647613260247372125243616270670;
const static double reb_whfast_corrector_b_115 = 0.072593394748842738674253180742744961827622366521517;

// Fast inverse factorial lookup table
static const double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

static inline double fastabs(double x){
        return (x > 0.) ? x : -x;
}

static void stumpff_cs(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fastabs(z)>0.1){
        z = z/4.;
        n++;
    }
    const int nmax = 15;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np=nmax-2;np>=5;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }
    cs[5] = c_odd;
    cs[4] = c_even;
    cs[3] = invfactorial[3]  - z *cs[5];
    cs[2] = invfactorial[2]  - z *cs[4];
    cs[1] = invfactorial[1]  - z *cs[3];
    for (;n>0;n--){ 
        z = z*4.;
        cs[5] = (cs[5]+cs[4]+cs[3]*cs[2])*0.0625;
        cs[4] = (1.+cs[1])*cs[3]*0.125;
        cs[3] = 1./6.-z*cs[5];
        cs[2] = 0.5-z*cs[4];
        cs[1] = 1.-z*cs[3];
    }
    cs[0] = invfactorial[0]  - z *cs[2];
}
static void stumpff_cs3(double *restrict cs, double z) {
    unsigned int n = 0;
    while(fabs(z)>0.1){
        z = z/4.;
        n++;
    }
    const int nmax = 13;
    double c_odd  = invfactorial[nmax];
    double c_even = invfactorial[nmax-1];
    for(int np=nmax-2;np>=3;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }
    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1]  - z *c_odd;
    cs[0] = invfactorial[0]  - z *c_even;
    for (;n>0;n--){ 
        cs[3] = (cs[2]+cs[0]*cs[3])*0.25;
        cs[2] = cs[1]*cs[1]*0.5;
        cs[1] = cs[0]*cs[1];
        cs[0] = 2.*cs[0]*cs[0]-1.;
    }
}

static void stiefel_Gs(double *restrict Gs, double beta, double X) {
    double X2 = X*X;
    stumpff_cs(Gs, beta*X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    double _pow = X2*X;
    Gs[3] *= _pow; 
    _pow *= X;
    Gs[4] *= _pow; 
    _pow *= X;
    Gs[5] *= _pow; 
    return;
}

static void stiefel_Gs3(double *restrict Gs, double beta, double X) {
    double X2 = X*X;
    stumpff_cs3(Gs, beta*X2);
    Gs[1] *= X; 
    Gs[2] *= X2; 
    Gs[3] *= X2*X;
    return;
}

#define WHFAST_NMAX_QUART 64    ///< Maximum number of iterations for quartic solver
#define WHFAST_NMAX_NEWT  32    ///< Maximum number of iterations for Newton's method
/************************************
 * Keplerian motion for one planet  */
void reb_whfast_kepler_solver(const struct reb_simulation* const r, struct reb_particle* const restrict p_j, const double M, unsigned int i, double _dt){
    const struct reb_particle p1 = p_j[i];

    const double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
    const double r0i = 1./r0;
    const double v2 =  p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
    const double beta = 2.*M*r0i - v2;
    const double eta0 = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;
    const double zeta0 = M - beta*r0;
    double X;
    double Gs[6]; 
    double invperiod=0;  // only used for beta>0. Set to 0 only to suppress compiler warnings.
    double X_per_period = nan(""); // only used for beta>0. nan triggers Newton's method for beta<0.
        
    if (beta>0.){
        // Elliptic orbit
        const double sqrt_beta = sqrt(beta);
        invperiod = sqrt_beta*beta/(2.*M_PI*M);
        X_per_period = 2.*M_PI/sqrt_beta;
        if (fabs(_dt)*invperiod>1. && r->ri_whfast.timestep_warning == 0){
            // Ignoring const qualifiers. This warning should not have any effect on
            // other parts of the code, nor is it vital to show it.
            ((struct reb_simulation* const)r)->ri_whfast.timestep_warning++;
            reb_warning((struct reb_simulation* const)r,"WHFast convergence issue. Timestep is larger than at least one orbital period.");
        }
        //X = _dt*invperiod*X_per_period; // first order guess 
        const double dtr0i = _dt*r0i;
        //X = dtr0i; // first order guess
        X = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess
        //X = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
        //X = _dt*beta/M + eta0/M*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version 
    }else{
        // Hyperbolic orbit
        X = 0.; // Initial guess 
    }

    unsigned int converged = 0;
    double oldX = X; 

    // Do one Newton step
    stiefel_Gs3(Gs, beta, X);
    const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
    double ri = 1./(r0 + eta0Gs1zeta0Gs2);
    X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);

    // Choose solver depending on estimated step size
    // Note, for hyperbolic orbits this uses Newton's method.
    if(fastabs(X-oldX) > 0.01*X_per_period){
        // Quartic solver
        // Linear initial guess
        X = beta*_dt/M;
        static double prevX[WHFAST_NMAX_QUART+1];
        for(int n_lag=1; n_lag < WHFAST_NMAX_QUART; n_lag++){
            stiefel_Gs3(Gs, beta, X);
            const double f = r0*X + eta0*Gs[2] + zeta0*Gs[3] - _dt;
            const double fp = r0 + eta0*Gs[1] + zeta0*Gs[2];
            const double fpp = eta0*Gs[0] + zeta0*Gs[1];
            const double denom = fp + sqrt(fabs(16.*fp*fp - 20.*f*fpp));
            X = (X*denom - 5.*f)/denom;
            for(int i=1;i<n_lag;i++){
                if(X==prevX[i]){
                    // Converged. Exit.
                    n_lag = WHFAST_NMAX_QUART;
                    converged = 1;
                    break;
                }
            }
            prevX[n_lag] = X;
        }
        const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
        ri = 1./(r0 + eta0Gs1zeta0Gs2);
    }else{
        // Newton's method
        double oldX2 = nan("");             
        for (int n_hg=1;n_hg<WHFAST_NMAX_NEWT;n_hg++){
            oldX2 = oldX;
            oldX = X;
            stiefel_Gs3(Gs, beta, X);
            const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
            ri = 1./(r0 + eta0Gs1zeta0Gs2);
            X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+_dt);
            
            if (X==oldX||X==oldX2){
                // Converged. Exit.
                converged = 1;
                break; 
            }
        }
    }
        
    // If solver did not work, fallback to bisection 
    if (converged == 0){ 
        double X_min, X_max;
        if (beta>0.){
            //Elliptic
            X_min = X_per_period * floor(_dt*invperiod);
            X_max = X_min + X_per_period;
        }else{
            //Hyperbolic
            double h2 = r0*r0*v2-eta0*eta0;
            double q = h2/M/(1.+sqrt(1.-h2*beta/(M*M)));
            double vq = copysign( sqrt(h2)/q, _dt);
            // X_max and X_min correspond to dt/r_min and dt/r_max
            // which are reachable in this timestep
            // r_max = vq*_dt+r0
            // r_min = pericenter
            X_min = _dt/(fastabs(vq*_dt)+r0); 
            X_max = _dt/q;
            if (_dt<0.){
                double temp = X_min;
                X_min = X_max;
                X_max = temp;
            }
        }
        X = (X_max + X_min)/2.;
        do{
            stiefel_Gs3(Gs, beta, X);
            double s   = r0*X + eta0*Gs[2] + zeta0*Gs[3]-_dt;
            if (s>=0.){
                X_max = X;
            }else{
                X_min = X;
            }
            X = (X_max + X_min)/2.;
        }while (fastabs((X_max-X_min))>fastabs((X_max+X_min)*1e-15));
        const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
        ri = 1./(r0 + eta0Gs1zeta0Gs2);
    }
    if (isnan(ri)){
        // Exception for (almost) straight line motion in hyperbolic case
        ri = 0.;
        Gs[1] = 0.;
        Gs[2] = 0.;
        Gs[3] = 0.;
    }

    // Note: These are not the traditional f and g functions.
    double f = -M*Gs[2]*r0i;
    double g = _dt - M*Gs[3];
    double fd = -M*Gs[1]*r0i*ri; 
    double gd = -M*Gs[2]*ri; 
        
    p_j[i].x += f*p1.x + g*p1.vx;
    p_j[i].y += f*p1.y + g*p1.vy;
    p_j[i].z += f*p1.z + g*p1.vz;
        
    p_j[i].vx += fd*p1.x + gd*p1.vx;
    p_j[i].vy += fd*p1.y + gd*p1.vy;
    p_j[i].vz += fd*p1.z + gd*p1.vz;

    //Variations
    for (int v=0;v<r->var_config_N;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        const int index = vc.index;
        stiefel_Gs(Gs, beta, X);    // Recalculate (to get Gs[4] and Gs[5])
        struct reb_particle dp1 = p_j[i+index];
        double dr0 = (dp1.x*p1.x + dp1.y*p1.y + dp1.z*p1.z)*r0i;
        double dbeta = -2.*M*dr0*r0i*r0i - 2.* (dp1.vx*p1.vx + dp1.vy*p1.vy + dp1.vz*p1.vz);
        double deta0 = dp1.x*p1.vx + dp1.y*p1.vy + dp1.z*p1.vz
                 + p1.x*dp1.vx + p1.y*dp1.vy + p1.z*dp1.vz;
        double dzeta0 = -beta*dr0 - r0*dbeta;
        double G3beta = 0.5*(3.*Gs[5]-X*Gs[4]);
        double G2beta = 0.5*(2.*Gs[4]-X*Gs[3]);
        double G1beta = 0.5*(Gs[3]-X*Gs[2]);
        double tbeta = eta0*G2beta + zeta0*G3beta;
        double dX = -1.*ri*(X*dr0 + Gs[2]*deta0+Gs[3]*dzeta0+tbeta*dbeta);
        double dG1 = Gs[0]*dX + G1beta*dbeta; 
        double dG2 = Gs[1]*dX + G2beta*dbeta;
        double dG3 = Gs[2]*dX + G3beta*dbeta;
        double dr = dr0 + Gs[1]*deta0 + Gs[2]*dzeta0 + eta0*dG1 + zeta0*dG2;
        double df = M*Gs[2]*dr0*r0i*r0i - M*dG2*r0i;
        double dg = -M*dG3;
        double dfd = -M*dG1*r0i*ri + M*Gs[1]*(dr0*r0i+dr*ri)*r0i*ri;
        double dgd = -M*dG2*ri + M*Gs[2]*dr*ri*ri;
    
        p_j[i+index].x += f*dp1.x + g*dp1.vx + df*p1.x + dg*p1.vx;
        p_j[i+index].y += f*dp1.y + g*dp1.vy + df*p1.y + dg*p1.vy;
        p_j[i+index].z += f*dp1.z + g*dp1.vz + df*p1.z + dg*p1.vz;

        p_j[i+index].vx += fd*dp1.x + gd*dp1.vx + dfd*p1.x + dgd*p1.vx;
        p_j[i+index].vy += fd*dp1.y + gd*dp1.vy + dfd*p1.y + dgd*p1.vy;
        p_j[i+index].vz += fd*dp1.z + gd*dp1.vz + dfd*p1.z + dgd*p1.vz;
    }

}

/***************************** 
 * Interaction Hamiltonian  */

void reb_whfast_interaction_step(struct reb_simulation* const r, const double _dt){
    const unsigned int N_real = r->N-r->N_var;
    const double G = r->G;
    struct reb_particle* particles = r->particles;
    const double m0 = particles[0].m;
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_j = ri_whfast->p_jh;
    switch (ri_whfast->coordinates){
        case REB_WHFAST_COORDINATES_JACOBI:
            {
            const double softening = r->softening;
            for (int v=0;v<r->var_config_N;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_inertial_to_jacobi_acc(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
            }
            reb_transformations_inertial_to_jacobi_acc(r->particles, p_j, r->particles, N_real);
            for (int v=0;v<r->var_config_N;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_inertial_to_jacobi_acc(r->particles+vc.index, p_j+vc.index, r->particles, N_real);
            }
            double eta = m0;
            for (unsigned int i=1;i<N_real;i++){
                // Eq 132
                const struct reb_particle pji = p_j[i];
                eta += pji.m;
                static double rj2i;
                static double rj3iM;
                static double prefac1;
                p_j[i].vx += _dt * pji.ax;
                p_j[i].vy += _dt * pji.ay;
                p_j[i].vz += _dt * pji.az;
                if (i>1){
                    rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z + softening*softening);
                    const double rji  = sqrt(rj2i);
                    rj3iM = rji*rj2i*G*eta;
                    prefac1 = _dt*rj3iM;
                    p_j[i].vx += prefac1*pji.x;
                    p_j[i].vy += prefac1*pji.y;
                    p_j[i].vz += prefac1*pji.z;
                }
                for(int v=0;v<r->var_config_N;v++){
                    struct reb_variational_configuration const vc = r->var_config[v];
                    const int index = vc.index;
                    double rj5M = rj3iM*rj2i;
                    double rdr = p_j[i+index].x*pji.x + p_j[i+index].y*pji.y + p_j[i+index].z*pji.z;
                    double prefac2 = -_dt*3.*rdr*rj5M;
                    p_j[i+index].vx += _dt * p_j[i+index].ax;
                    p_j[i+index].vy += _dt * p_j[i+index].ay;
                    p_j[i+index].vz += _dt * p_j[i+index].az;
                    if (i>1){
                        p_j[i+index].vx += prefac1*p_j[i+index].x + prefac2*pji.x;
                        p_j[i+index].vy += prefac1*p_j[i+index].y + prefac2*pji.y;
                        p_j[i+index].vz += prefac1*p_j[i+index].z + prefac2*pji.z;
                    }
                }
            }
            }
            break;
        case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
#pragma omp parallel for 
            for (unsigned int i=1;i<N_real;i++){
                p_j[i].vx += _dt*particles[i].ax;
                p_j[i].vy += _dt*particles[i].ay;
                p_j[i].vz += _dt*particles[i].az;
            }
            break;
        case REB_WHFAST_COORDINATES_WHDS:
#pragma omp parallel for 
            for (unsigned int i=1;i<N_real;i++){
                const double mi = particles[i].m;
                p_j[i].vx += _dt*(m0+mi)*particles[i].ax/m0;
                p_j[i].vy += _dt*(m0+mi)*particles[i].ay/m0;
                p_j[i].vz += _dt*(m0+mi)*particles[i].az/m0;
            }
            break;
    };
}
void reb_whfast_jump_step(const struct reb_simulation* const r, const double _dt){
    const struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* const p_h = r->ri_whfast.p_jh;
    const int N_real = r->N-r->N_var;
    const double m0 = r->particles[0].m;
    switch (ri_whfast->coordinates){
        case REB_WHFAST_COORDINATES_JACOBI:
            // Nothing to be done.
            break;
        case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            {
            double px=0, py=0, pz=0;
#pragma omp parallel for reduction (+:px), reduction (+:py), reduction (+:pz)
            for(int i=1;i<N_real;i++){
                const double m = r->particles[i].m;
                px += m * p_h[i].vx;
                py += m * p_h[i].vy;
                pz += m * p_h[i].vz;
            }
#pragma omp parallel for 
            for(int i=1;i<N_real;i++){
                p_h[i].x += _dt * (px/m0);
                p_h[i].y += _dt * (py/m0);
                p_h[i].z += _dt * (pz/m0);
            }
            }
            break;
        case REB_WHFAST_COORDINATES_WHDS:
            {
            double px=0, py=0, pz=0;
#pragma omp parallel for reduction (+:px), reduction (+:py), reduction (+:pz)
            for(int i=1;i<N_real;i++){
                const double m = r->particles[i].m;
                px += m * p_h[i].vx / (m0+m);
                py += m * p_h[i].vy / (m0+m);
                pz += m * p_h[i].vz / (m0+m);
             }
#pragma omp parallel for 
            for(int i=1;i<N_real;i++){
                const double m = r->particles[i].m;
                p_h[i].x += _dt * (px - (m * p_h[i].vx / (m0+m)) );
                p_h[i].y += _dt * (py - (m * p_h[i].vy / (m0+m)) );
                p_h[i].z += _dt * (pz - (m * p_h[i].vz / (m0+m)) );
            }
            }
            break;
    };
}

/***************************** 
 * DKD Scheme                */

void reb_whfast_kepler_step(const struct reb_simulation* const r, const double _dt){
    const double m0 = r->particles[0].m;
    const double G = r->G;
    const unsigned int N_real = r->N-r->N_var;
    const int coordinates = r->ri_whfast.coordinates;
    struct reb_particle* const p_j = r->ri_whfast.p_jh;
    double eta = m0;
#pragma omp parallel for 
    for (unsigned int i=1;i<N_real;i++){
        switch (coordinates){
            case REB_WHFAST_COORDINATES_JACOBI:
                eta += p_j[i].m;
                break;
            case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                //  eta = m0
                break;
            case REB_WHFAST_COORDINATES_WHDS:
                eta = m0+p_j[i].m;
                break;
        };
        reb_whfast_kepler_solver(r, p_j, eta*G, i, _dt);
    }
}

void reb_whfast_com_step(const struct reb_simulation* const r, const double _dt){
    struct reb_particle* const p_j = r->ri_whfast.p_jh;
    p_j[0].x += _dt*p_j[0].vx;
    p_j[0].y += _dt*p_j[0].vy;
    p_j[0].z += _dt*p_j[0].vz;
}

static void reb_whfast_corrector_Z(struct reb_simulation* r, const double a, const double b){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N-r->N_var;
    reb_whfast_kepler_step(r, a);
    reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
    for (int v=0;v<r->var_config_N;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        reb_transformations_jacobi_to_inertial_pos(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
    }
    r->gravity_ignore_terms = 1;
    reb_update_acceleration(r);
    reb_whfast_interaction_step(r, -b);
    reb_whfast_kepler_step(r, -2.*a);
    reb_transformations_jacobi_to_inertial_pos(particles, ri_whfast->p_jh, particles, N_real);
    for (int v=0;v<r->var_config_N;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        reb_transformations_jacobi_to_inertial_pos(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
    }
    r->gravity_ignore_terms = 1;
    reb_update_acceleration(r);
    reb_whfast_interaction_step(r, b);
    reb_whfast_kepler_step(r, a);
}

void reb_whfast_apply_corrector(struct reb_simulation* r, double inv, int order, void (*corrector_Z)(struct reb_simulation*, const double, const double)){
    const double dt = r->dt;
    if (order==3){
        // Third order corrector
        corrector_Z(r, reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_31*dt);
        corrector_Z(r, -reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_31*dt);
    }
    if (order==5){
        // Fifth order corrector
        corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_51*dt);
        corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_52*dt);
        corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_52*dt);
        corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_51*dt);
    }
    if (order==7){
        // Seventh order corrector
        corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_71*dt);
        corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_72*dt);
        corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_73*dt);
        corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_73*dt);
        corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_72*dt);
        corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_71*dt);
    }
    if (order==11){
        // Eleventh order corrector
        corrector_Z(r, -reb_whfast_corrector_a_5*dt,-inv*reb_whfast_corrector_b_111*dt);
        corrector_Z(r, -reb_whfast_corrector_a_4*dt,-inv*reb_whfast_corrector_b_112*dt);
        corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_113*dt);
        corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_114*dt);
        corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_115*dt);
        corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_115*dt);
        corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_114*dt);
        corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_113*dt);
        corrector_Z(r, reb_whfast_corrector_a_4*dt,inv*reb_whfast_corrector_b_112*dt);
        corrector_Z(r, reb_whfast_corrector_a_5*dt,inv*reb_whfast_corrector_b_111*dt);
    }
}

int reb_integrator_whfast_init(struct reb_simulation* const r){
    for (int v=0;v<r->var_config_N;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        if (vc.order!=1){
            reb_error(r, "WHFast/MEGNO only supports first order variational equations.");
            return 1; // Error
        }
        if (vc.testparticle>=0){
            reb_error(r, "Test particle variations not supported with WHFast. Use IAS15.");
            return 1; // Error
        }
    }
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
#if defined(_OPENMP)
    if (ri_whfast->coordinates!=REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC
        && ri_whfast->coordinates!=REB_WHFAST_COORDINATES_WHDS){
        reb_error(r,"WHFast when used with OpenMP requires REB_WHFAST_COORDINATES_WHDS or REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC\n");
        return 1; // Error
    }
#endif
    if (r->var_config_N>0 && ri_whfast->coordinates!=REB_WHFAST_COORDINATES_JACOBI){
        reb_error(r, "Variational particles are only compatible with Jacobi coordinates.");
        return 1; // Error
    }
    if (ri_whfast->corrector!=0 && ri_whfast->coordinates!=REB_WHFAST_COORDINATES_JACOBI){
        reb_error(r, "Symplectic correctors are only compatible with Jacobi coordinates.");
        return 1; // Error
    }
    const int N = r->N;
    if (ri_whfast->coordinates==REB_WHFAST_COORDINATES_JACOBI){
        r->gravity_ignore_terms = 1;
    }else{
        r->gravity_ignore_terms = 2;
    }
    if (ri_whfast->allocated_N != N){
        ri_whfast->allocated_N = N;
        ri_whfast->p_jh = realloc(ri_whfast->p_jh,sizeof(struct reb_particle)*N);
        ri_whfast->recalculate_coordinates_this_timestep = 1;
    }
    return 0;
}

void reb_integrator_whfast_from_inertial(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    const int N_real = N-r->N_var;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    
    switch (ri_whfast->coordinates){
        case REB_WHFAST_COORDINATES_JACOBI:
            reb_transformations_inertial_to_jacobi_posvel(particles, ri_whfast->p_jh, particles, N_real);
            for (int v=0;v<r->var_config_N;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_inertial_to_jacobi_posvel(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
            }
            break;
        case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            reb_transformations_inertial_to_democraticheliocentric_posvel_testparticles(particles, ri_whfast->p_jh, N_real, N_active);
            break;
        case REB_WHFAST_COORDINATES_WHDS:
            reb_transformations_inertial_to_whds_posvel(particles, ri_whfast->p_jh, N_real);
            break;
    };
}

void reb_integrator_whfast_to_inertial(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    const int N_real = N-r->N_var;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    
    // Prepare coordinates for KICK step
    if (r->force_is_velocity_dependent){
        switch (ri_whfast->coordinates){
            case REB_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(particles, ri_whfast->p_jh, particles, N_real);
                break;
            case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel_testparticles(particles, ri_whfast->p_jh, N_real, N_active);
                break;
            case REB_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(particles, ri_whfast->p_jh, N_real);
                break;
        };
    }else{
        switch (ri_whfast->coordinates){
            case REB_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(particles, ri_whfast->p_jh, particles, N_real);
                break;
            case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel_testparticles(particles, ri_whfast->p_jh, N_real, N_active);
                break;
            case REB_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(particles, ri_whfast->p_jh, N_real);
                break;
        };
    }
}


void reb_integrator_whfast_part1(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    const int N_real = N-r->N_var;
    
    if (reb_integrator_whfast_init(r)){
        // Non recoverable error occured.
        return;
    }
    
    // Only recalculate Jacobi coordinates if needed
    if (ri_whfast->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
        if (ri_whfast->is_synchronized==0){
            reb_integrator_whfast_synchronize(r);
            if (ri_whfast->recalculate_coordinates_but_not_synchronized_warning==0){
                reb_warning(r,"Recalculating coordinates but pos/vel were not synchronized before.");
                ri_whfast->recalculate_coordinates_but_not_synchronized_warning++;
            }
        }
        reb_integrator_whfast_from_inertial(r);
        ri_whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (ri_whfast->is_synchronized){
        // First half DRIFT step
        if (ri_whfast->corrector){
            reb_whfast_apply_corrector(r, 1., ri_whfast->corrector, reb_whfast_corrector_Z);
        }
        reb_whfast_kepler_step(r, r->dt/2.);    // half timestep
        reb_whfast_com_step(r, r->dt/2.);
    }else{
        // Combined DRIFT step
        reb_whfast_kepler_step(r, r->dt);    // full timestep
        reb_whfast_com_step(r, r->dt);
    }

    reb_whfast_jump_step(r,r->dt/2.);

    reb_integrator_whfast_to_inertial(r);

    // Variational equations only available for jacobi coordinates. 
    // If other coordinates are used, the code will raise an exception in part1 of the integrator.
    for (int v=0;v<r->var_config_N;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        ri_whfast->p_jh[vc.index].x += r->dt/2.*ri_whfast->p_jh[vc.index].vx;
        ri_whfast->p_jh[vc.index].y += r->dt/2.*ri_whfast->p_jh[vc.index].vy;
        ri_whfast->p_jh[vc.index].z += r->dt/2.*ri_whfast->p_jh[vc.index].vz;
        if (r->force_is_velocity_dependent){
            reb_transformations_jacobi_to_inertial_posvel(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
        }else{
            reb_transformations_jacobi_to_inertial_pos(particles+vc.index, ri_whfast->p_jh+vc.index, particles, N_real);
        }
    }

    r->t+=r->dt/2.;
}

void reb_integrator_whfast_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    if (ri_whfast->is_synchronized == 0){
        const int N_real = r->N-r->N_var;
        const int N_active = r->N_active==-1?r->N:r->N_active;
        struct reb_particle* sync_pj  = NULL;
        if (ri_whfast->keep_unsynchronized){
            sync_pj = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(sync_pj,r->ri_whfast.p_jh,r->N*sizeof(struct reb_particle));
        }
        reb_whfast_kepler_step(r, r->dt/2.);
        reb_whfast_com_step(r, r->dt/2.);
        if (ri_whfast->corrector){
            reb_whfast_apply_corrector(r, -1., ri_whfast->corrector, reb_whfast_corrector_Z);
        }
        switch (ri_whfast->coordinates){
            case REB_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(r->particles, ri_whfast->p_jh, r->particles, N_real);
                break;
            case REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel_testparticles(r->particles, ri_whfast->p_jh, N_real, N_active);
                break;
            case REB_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(r->particles, ri_whfast->p_jh, N_real);
                break;
        };
        for (int v=0;v<r->var_config_N;v++){
            struct reb_variational_configuration const vc = r->var_config[v];
            reb_transformations_jacobi_to_inertial_posvel(r->particles+vc.index, ri_whfast->p_jh+vc.index, r-> particles, N_real);
        }
        if (ri_whfast->keep_unsynchronized){
            memcpy(r->ri_whfast.p_jh,sync_pj,r->N*sizeof(struct reb_particle));
            free(sync_pj);
        }else{
            ri_whfast->is_synchronized = 1;
        }
    }
}

void reb_integrator_whfast_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    if (ri_whfast->p_jh==NULL){
        // Non recoverable error occured earlier. 
        // Skipping rest of integration to avoid segmentation fault.
        return;
    }
    
    reb_whfast_interaction_step(r, r->dt);
    reb_whfast_jump_step(r,r->dt/2.);
   
    ri_whfast->is_synchronized = 0;
    if (ri_whfast->safe_mode){
        reb_integrator_whfast_synchronize(r);
    }
    
    r->t+=r->dt/2.;
    r->dt_last_done = r->dt;

    
    if (r->var_config_N){
        // Need to have x,v,a synchronized to calculate ddot/d for MEGNO. 
        reb_integrator_whfast_synchronize(r);
        // Add additional acceleration term for MEGNO calculation
        const int N_real = r->N-r->N_var;
        struct reb_particle* restrict const particles = r->particles;
        for (int v=0;v<r->var_config_N;v++){
            struct reb_variational_configuration const vc = r->var_config[v];
            struct reb_particle* const particles_var1 = particles + vc.index;
            const int index = vc.index;
            // Centre of mass
            ri_whfast->p_jh[index].x += r->dt/2.*ri_whfast->p_jh[index].vx;
            ri_whfast->p_jh[index].y += r->dt/2.*ri_whfast->p_jh[index].vy;
            ri_whfast->p_jh[index].z += r->dt/2.*ri_whfast->p_jh[index].vz;
            reb_transformations_jacobi_to_inertial_posvel(particles_var1, ri_whfast->p_jh+index, particles, N_real);
            if (r->calculate_megno){
                reb_calculate_acceleration_var(r);
                const double dx = particles[0].x - particles[1].x;
                const double dy = particles[0].y - particles[1].y;
                const double dz = particles[0].z - particles[1].z;
                const double r2 = dx*dx + dy*dy + dz*dz + r->softening*r->softening;
                const double _r  = sqrt(r2);
                const double r3inv = 1./(r2*_r);
                const double r5inv = 3.*r3inv/r2;
                const double ddx = particles_var1[0].x - particles_var1[1].x;
                const double ddy = particles_var1[0].y - particles_var1[1].y;
                const double ddz = particles_var1[0].z - particles_var1[1].z;
                const double Gmi = r->G * particles[0].m;
                const double Gmj = r->G * particles[1].m;
                const double dax =   ddx * ( dx*dx*r5inv - r3inv )
                           + ddy * ( dx*dy*r5inv )
                           + ddz * ( dx*dz*r5inv );
                const double day =   ddx * ( dy*dx*r5inv )
                           + ddy * ( dy*dy*r5inv - r3inv )
                           + ddz * ( dy*dz*r5inv );
                const double daz =   ddx * ( dz*dx*r5inv )
                           + ddy * ( dz*dy*r5inv )
                           + ddz * ( dz*dz*r5inv - r3inv );
                
                particles_var1[0].ax += Gmj * dax;
                particles_var1[0].ay += Gmj * day;
                particles_var1[0].az += Gmj * daz;
                
                particles_var1[1].ax -= Gmi * dax;
                particles_var1[1].ay -= Gmi * day;
                particles_var1[1].az -= Gmi * daz;

                // TODO Need to add mass terms. Also need to add them to tangent map above.
            }
        }

        // Update MEGNO in middle of timestep as we need synchonized x/v/a.
        if (r->calculate_megno){
            double dY = r->dt * 2. * r->t * reb_tools_megno_deltad_delta(r);
            reb_tools_megno_update(r, dY);
        }
    }
}
    
void reb_integrator_whfast_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfast* const ri_whfast = &(r->ri_whfast);
    ri_whfast->corrector = 0;
    ri_whfast->coordinates = REB_WHFAST_COORDINATES_JACOBI;
    ri_whfast->is_synchronized = 1;
    ri_whfast->keep_unsynchronized = 0;
    ri_whfast->safe_mode = 1;
    ri_whfast->recalculate_coordinates_this_timestep = 0;
    ri_whfast->allocated_N = 0;
    ri_whfast->timestep_warning = 0;
    ri_whfast->recalculate_coordinates_but_not_synchronized_warning = 0;
    if (ri_whfast->p_jh){
        free(ri_whfast->p_jh);
        ri_whfast->p_jh = NULL;
    }
}
