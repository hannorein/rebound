/**
 * @file    integrator_whfast.c
 * @brief   WHFAST integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the WHFast integration scheme.  
 * Described in Rein & Tamayo 2015, Rein et al. 2019.
 * See also Wisdom et al. 1996.
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

#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include <math.h>
#include <stddef.h>
#include "transformations.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "binarydata.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

const struct reb_binarydata_field_descriptor reb_integrator_whfast_field_descriptor_list[] = {
    { 61, REB_UINT,         "corrector",          offsetof(struct reb_integrator_whfast_state, corrector), 0, 0, 0},
    { 62, REB_UINT,         "recalculate_coordinates_this_timestep", offsetof(struct reb_integrator_whfast_state, recalculate_coordinates_this_timestep), 0, 0, 0},
    { 63, REB_UINT,         "safe_mode",          offsetof(struct reb_integrator_whfast_state, safe_mode), 0, 0, 0},
    { 64, REB_UINT,         "keep_unsynchronized",offsetof(struct reb_integrator_whfast_state, keep_unsynchronized), 0, 0, 0},
    { 104, REB_POINTER,     "p_jh",               offsetof(struct reb_integrator_whfast_state, p_jh), offsetof(struct reb_integrator_whfast_state, N_allocated), sizeof(struct reb_particle), 0},
    { 105, REB_POINTER,     "p_jh_var",           offsetof(struct reb_integrator_whfast_state, p_jh_var), offsetof(struct reb_integrator_whfast_state, N_allocated_var), sizeof(struct reb_particle), 0},
    { 117, REB_INT,         "coordinates",        offsetof(struct reb_integrator_whfast_state, coordinates), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_WHFAST_COORDINATES) },
    { 143, REB_UINT,        "corrector2",         offsetof(struct reb_integrator_whfast_state, corrector2), 0, 0, 0},
    { 144, REB_INT,         "kernel",             offsetof(struct reb_integrator_whfast_state, kernel), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_WHFAST_KERNEL)},
    { 0 }, // Null terminated list
};

void* reb_integrator_whfast_create();
void reb_integrator_whfast_free(void* state);
void reb_integrator_whfast_synchronize(struct reb_simulation* const r, void* state);
void reb_integrator_whfast_step(struct reb_simulation* const r, void* state);

const struct reb_integrator reb_integrator_whfast = {
    .step = reb_integrator_whfast_step,
    .synchronize = reb_integrator_whfast_synchronize,
    .create = reb_integrator_whfast_create,
    .free = reb_integrator_whfast_free,
    .field_descriptor_list = reb_integrator_whfast_field_descriptor_list,
};

void* reb_integrator_whfast_create(){
    struct reb_integrator_whfast_state* whfast = calloc(sizeof(struct reb_integrator_whfast_state),1);
    whfast->corrector = 0;
    whfast->corrector2 = 0;
    whfast->kernel = 0;
    whfast->coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI;
    whfast->keep_unsynchronized = 0;
    whfast->safe_mode = 1;
    whfast->recalculate_coordinates_this_timestep = 0;
    whfast->recalculate_coordinates_but_not_synchronized_warning = 0;
    return whfast;
}

void reb_integrator_whfast_free(void* p){
    struct reb_integrator_whfast_state* whfast = p;
    free(whfast->p_jh);
    free(whfast->p_jh_var);
    free(whfast->p_temp);
    free(whfast);
}

// Corrector coefficients
static const double reb_whfast_corrector_a_1 = 0.41833001326703777398908601289259374469640768464934;
static const double reb_whfast_corrector_a_2 = 0.83666002653407554797817202578518748939281536929867;
static const double reb_whfast_corrector_a_3 = 1.2549900398011133219672580386777812340892230539480;
static const double reb_whfast_corrector_a_4 = 1.6733200530681510959563440515703749787856307385973;
static const double reb_whfast_corrector_a_5 = 2.0916500663351888699454300644629687234820384232467;
static const double reb_whfast_corrector_a_6 = 2.5099800796022266439345160773555624681784461078960; 
static const double reb_whfast_corrector_a_7 = 2.9283100928692644179236020902481562128748537925454;
static const double reb_whfast_corrector_a_8 = 3.3466401061363021919126881031407499575712614771947;
static const double reb_whfast_corrector_b_31 = -0.024900596027799867499350357910273437184309981229127;
static const double reb_whfast_corrector_b_51 = -0.0083001986759332891664501193034244790614366604097090;
static const double reb_whfast_corrector_b_52 = 0.041500993379666445832250596517122395307183302048545;
static const double reb_whfast_corrector_b_71 = 0.0024926811426922105779030593952776964450539008582219;
static const double reb_whfast_corrector_b_72 = -0.018270923246702131478062356884535264841652263842597;
static const double reb_whfast_corrector_b_73 = 0.053964399093127498721765893493510877532452806339655;
static const double reb_whfast_corrector_b_111 = 0.00020361579647854651301632818774633716473696537436847;
static const double reb_whfast_corrector_b_112 = -0.0023487215292295354188307328851055489876255097419754;
static const double reb_whfast_corrector_b_113 = 0.012309078592019946317544564763237909911330686448336;
static const double reb_whfast_corrector_b_114 = -0.038121613681288650508647613260247372125243616270670;
static const double reb_whfast_corrector_b_115 = 0.072593394748842738674253180742744961827622366521517;
static const double reb_whfast_corrector_b_178 = 0.093056103771425958591541059067553547100903397724386; 
static const double reb_whfast_corrector_b_177 = -0.065192863576377893658290760803725762027864651086787; 
static const double reb_whfast_corrector_b_176 = 0.032422198864713580293681523029577130832258806467604; 
static const double reb_whfast_corrector_b_175 = -0.012071760822342291062449751726959664253913904872527; 
static const double reb_whfast_corrector_b_174 = 0.0033132577069380655655490196833451994080066801611459; 
static const double reb_whfast_corrector_b_173 = -0.00063599983075817658983166881625078545864140848560259; 
static const double reb_whfast_corrector_b_172 = 0.000076436355227935738363241846979413475106795392377415; 
static const double reb_whfast_corrector_b_171 = -0.0000043347415473373580190650223498124944896789841432241; 
static const double reb_whfast_corrector2_b = 0.03486083443891981449909050107438281205803;

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
                                // Keplerian motion for one planet                       
                                // r only needed for variational particles and warning. Can be NULL.
void reb_integrator_whfast_kepler_solver(struct reb_particle* const restrict p, double mu, double dt, const struct reb_simulation* const r){
    const struct reb_particle p1 = *p; // Copy of particle

    const double r0 = sqrt(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z);
    const double r0i = 1./r0;
    const double v2 =  p1.vx*p1.vx + p1.vy*p1.vy + p1.vz*p1.vz;
    const double beta = 2.*mu*r0i - v2;
    const double eta0 = p1.x*p1.vx + p1.y*p1.vy + p1.z*p1.vz;
    const double zeta0 = mu - beta*r0;
    double X;
    double Gs[6]; 
    double invperiod=0;  // only used for beta>0. Set to 0 only to suppress compiler warnings.
    double X_per_period = nan(""); // only used for beta>0. nan triggers Newton's method for beta<0.

    if (beta>0.){
        // Elliptic orbit
        const double sqrt_beta = sqrt(beta);
        invperiod = sqrt_beta*beta/(2.*M_PI*mu);
        X_per_period = 2.*M_PI/sqrt_beta;
        if (fabs(dt)*invperiod>1.){
            if (r && !(r->messages_timestep_warning & 1)){
                ((struct reb_simulation*)r)->messages_timestep_warning |= 1; // casting away const qualifier is ok here.
                reb_simulation_warning(((struct reb_simulation*)r), "Possible convergence issue. Timestep in Kepler solver is larger than one orbital period.");
            }
        }

        //X = dt*invperiod*X_per_period; // first order guess 
        const double dtr0i = dt*r0i;
        //X = dtr0i; // first order guess
        X = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess
                                               //X = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
                                               //X = dt*beta/mu + eta0/mu*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version 
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
    X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+dt);

    // Choose solver depending on estimated step size
    // Note, for hyperbolic orbits this uses Newton's method.
    if(fastabs(X-oldX) > 0.01*X_per_period){
        // Quartic solver
        // Linear initial guess
        X = beta*dt/mu;
        double prevX[WHFAST_NMAX_QUART+1];
        for(int n_lag=1; n_lag < WHFAST_NMAX_QUART; n_lag++){
            stiefel_Gs3(Gs, beta, X);
            const double f = r0*X + eta0*Gs[2] + zeta0*Gs[3] - dt;
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
            X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+dt);

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
            X_min = X_per_period * floor(dt*invperiod);
            X_max = X_min + X_per_period;
        }else{
            //Hyperbolic
            double h2 = r0*r0*v2-eta0*eta0;
            double q = h2/mu/(1.+sqrt(1.-h2*beta/(mu*mu)));
            double vq = copysign( sqrt(h2)/q, dt);
            // X_max and X_min correspond to dt/r_min and dt/r_max
            // which are reachable in this timestep
            // r_max = vq*dt+r0
            // r_min = pericenter
            X_min = dt/(fastabs(vq*dt)+r0); 
            X_max = dt/q;
            if (dt<0.){
                double temp = X_min;
                X_min = X_max;
                X_max = temp;
            }
        }
        X = (X_max + X_min)/2.;
        do{
            stiefel_Gs3(Gs, beta, X);
            double s   = r0*X + eta0*Gs[2] + zeta0*Gs[3]-dt;
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
    double f = -mu*Gs[2]*r0i;
    double g = dt - mu*Gs[3];
    double fd = -mu*Gs[1]*r0i*ri; 
    double gd = -mu*Gs[2]*ri; 

    p->x += f*p1.x + g*p1.vx;
    p->y += f*p1.y + g*p1.vy;
    p->z += f*p1.z + g*p1.vz;

    p->vx += fd*p1.x + gd*p1.vx;
    p->vy += fd*p1.y + gd*p1.vy;
    p->vz += fd*p1.z + gd*p1.vz;

    //Variations
    for (size_t v=0;r && v<r->N_var_config;v++){
        struct reb_variational_configuration const vc = r->var_config[v];
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        const size_t pindex = p - whfast->p_jh;
        struct reb_particle* dp1p = whfast->p_jh_var + pindex + vc.index;
        struct reb_particle dp1 = *dp1p;
        stiefel_Gs(Gs, beta, X);    // Recalculate (to get Gs[4] and Gs[5])
        double dr0 = (dp1.x*p1.x + dp1.y*p1.y + dp1.z*p1.z)*r0i;
        double dbeta = -2.*mu*dr0*r0i*r0i - 2.* (dp1.vx*p1.vx + dp1.vy*p1.vy + dp1.vz*p1.vz);
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
        double df = mu*Gs[2]*dr0*r0i*r0i - mu*dG2*r0i;
        double dg = -mu*dG3;
        double dfd = -mu*dG1*r0i*ri + mu*Gs[1]*(dr0*r0i+dr*ri)*r0i*ri;
        double dgd = -mu*dG2*ri + mu*Gs[2]*dr*ri*ri;

        dp1p->x += f*dp1.x + g*dp1.vx + df*p1.x + dg*p1.vx;
        dp1p->y += f*dp1.y + g*dp1.vy + df*p1.y + dg*p1.vy;
        dp1p->z += f*dp1.z + g*dp1.vz + df*p1.z + dg*p1.vz;

        dp1p->vx += fd*dp1.x + gd*dp1.vx + dfd*p1.x + dgd*p1.vx;
        dp1p->vy += fd*dp1.y + gd*dp1.vy + dfd*p1.y + dgd*p1.vy;
        dp1p->vz += fd*dp1.z + gd*dp1.vz + dfd*p1.z + dgd*p1.vz;
    }
}

/***************************** 
 * Interaction Hamiltonian  */
void reb_integrator_whfast_interaction_step(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt){
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    const double G = r->G;
    struct reb_particle* particles = r->particles;
    struct reb_particle* particles_var = r->particles_var;
    const double m0 = particles[0].m;
    switch (coordinates){
        case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
            {
                reb_transformations_inertial_to_jacobi_acc(particles, p_jh, particles, N, N_active);
                for (size_t v=0;v<r->N_var_config;v++){
                    // Only WHFast supports variational equations
                    struct reb_integrator_whfast_state* whfast = r->integrator.state;
                    struct reb_variational_configuration const vc = r->var_config[v];
                    reb_transformations_inertial_to_jacobi_acc(particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
                }
                double eta = m0;
                for (size_t i=1;i<N;i++){
                    // Eq 132
                    const struct reb_particle pji = p_jh[i];
                    if (i<N_active){
                        eta += pji.m;
                    }
                    p_jh[i].vx += _dt * pji.ax;
                    p_jh[i].vy += _dt * pji.ay;
                    p_jh[i].vz += _dt * pji.az;
                    if (r->gravity != REB_GRAVITY_JACOBI){ 
                        // If Jacobi terms have not been added in update_acceleration, then add them here:
                        if (i>1){
                            const double rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
                            const double rji  = sqrt(rj2i);
                            const double rj3iM = rji*rj2i*G*eta;
                            const double prefac1 = _dt*rj3iM;
                            p_jh[i].vx += prefac1*pji.x;
                            p_jh[i].vy += prefac1*pji.y;
                            p_jh[i].vz += prefac1*pji.z;
                            for(size_t v=0;v<r->N_var_config;v++){
                                struct reb_integrator_whfast_state* whfast = r->integrator.state;
                                struct reb_particle* const p_jh_var = whfast->p_jh_var;
                                struct reb_variational_configuration const vc = r->var_config[v];
                                const size_t index = vc.index;
                                double rj5M = rj3iM*rj2i;
                                double rdr = p_jh_var[i+index].x*pji.x + p_jh_var[i+index].y*pji.y + p_jh_var[i+index].z*pji.z;
                                double prefac2 = -_dt*3.*rdr*rj5M;
                                p_jh_var[i+index].vx += prefac1*p_jh_var[i+index].x + prefac2*pji.x;
                                p_jh_var[i+index].vy += prefac1*p_jh_var[i+index].y + prefac2*pji.y;
                                p_jh_var[i+index].vz += prefac1*p_jh_var[i+index].z + prefac2*pji.z;
                            }
                        }
                        for(size_t v=0;v<r->N_var_config;v++){
                            struct reb_integrator_whfast_state* whfast = r->integrator.state;
                            struct reb_particle* const p_jh_var = whfast->p_jh_var;
                            struct reb_variational_configuration const vc = r->var_config[v];
                            const size_t index = vc.index;
                            p_jh_var[i+index].vx += _dt * p_jh_var[i+index].ax;
                            p_jh_var[i+index].vy += _dt * p_jh_var[i+index].ay;
                            p_jh_var[i+index].vz += _dt * p_jh_var[i+index].az;
                        }
                    }
                }
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
#pragma omp parallel for 
            for (size_t i=1;i<N;i++){
                p_jh[i].vx += _dt*particles[i].ax;
                p_jh[i].vy += _dt*particles[i].ay;
                p_jh[i].vz += _dt*particles[i].az;
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
#pragma omp parallel for 
            for (size_t i=1;i<N_active;i++){
                const double mi = particles[i].m;
                p_jh[i].vx += _dt*(m0+mi)*particles[i].ax/m0;
                p_jh[i].vy += _dt*(m0+mi)*particles[i].ay/m0;
                p_jh[i].vz += _dt*(m0+mi)*particles[i].az/m0;
            }
#pragma omp parallel for 
            for (size_t i=N_active;i<N;i++){
                p_jh[i].vx += _dt*particles[i].ax;
                p_jh[i].vy += _dt*particles[i].ay;
                p_jh[i].vz += _dt*particles[i].az;
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
            for (size_t i=1;i<N;i++){
                const double dr = sqrt(p_jh[i].x*p_jh[i].x + p_jh[i].y*p_jh[i].y + p_jh[i].z*p_jh[i].z);
                const double prefac = G*p_jh[0].m/(dr*dr*dr);
                p_jh[i].vx += _dt*(prefac*p_jh[i].x + particles[i].ax);
                p_jh[i].vy += _dt*(prefac*p_jh[i].y + particles[i].ay);
                p_jh[i].vz += _dt*(prefac*p_jh[i].z + particles[i].az);
            }
            break;
    };
}
void reb_integrator_whfast_jump_step(const struct reb_simulation* const r, struct reb_integrator_whfast_state* whfast, const double _dt){
    struct reb_particle* const p_h = whfast->p_jh;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    const double m0 = r->particles[0].m;
    switch (whfast->coordinates){
        case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
            // Nothing to be done.
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            {
                double px=0, py=0, pz=0;
#pragma omp parallel for reduction (+:px), reduction (+:py), reduction (+:pz)
                for(size_t i=1;i<N_active;i++){
                    const double m = r->particles[i].m;
                    px += m * p_h[i].vx;
                    py += m * p_h[i].vy;
                    pz += m * p_h[i].vz;
                }
#pragma omp parallel for 
                for(size_t i=1;i<N;i++){
                    p_h[i].x += _dt * (px/m0);
                    p_h[i].y += _dt * (py/m0);
                    p_h[i].z += _dt * (pz/m0);
                }
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
            {
                double px=0, py=0, pz=0;
#pragma omp parallel for reduction (+:px), reduction (+:py), reduction (+:pz)
                for(size_t i=1;i<N_active;i++){
                    const double m = r->particles[i].m;
                    px += m * p_h[i].vx / (m0+m);
                    py += m * p_h[i].vy / (m0+m);
                    pz += m * p_h[i].vz / (m0+m);
                }
#pragma omp parallel for 
                for(size_t i=1;i<N_active;i++){
                    const double m = r->particles[i].m;
                    p_h[i].x += _dt * (px - (m * p_h[i].vx / (m0+m)) );
                    p_h[i].y += _dt * (py - (m * p_h[i].vy / (m0+m)) );
                    p_h[i].z += _dt * (pz - (m * p_h[i].vz / (m0+m)) );
                }
#pragma omp parallel for 
                for(size_t i=N_active;i<N;i++){
                    p_h[i].x += _dt * px;
                    p_h[i].y += _dt * py;
                    p_h[i].z += _dt * pz;
                }
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
            // Nothing to be done.
            break;
    };
}

/***************************** 
 * DKD Scheme                */

void reb_integrator_whfast_kepler_step(const struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt){
    const double m0 = r->particles[0].m;
    const double G = r->G;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    double eta = m0;
    switch (coordinates){
        case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
#pragma omp parallel for private(eta)
            for (size_t i=1;i<N;i++){
#ifdef OPENMP
                eta = m0;
                for (size_t j=1;j<MIN(i,N_active);j++){
                    eta += p_jh[j].m;
                }
#else // OPENMP
                if (i<N_active){
                    eta += p_jh[i].m;
                }
#endif // OPENMP
                reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r);
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
#pragma omp parallel for 
            for (size_t i=1;i<N;i++){
                reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r); // eta = m0
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
#pragma omp parallel for private(eta)
            for (size_t i=1;i<N;i++){
                if (i<N_active){
                    eta = m0+p_jh[i].m;
                }else{
                    eta = m0;
                }
                reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r);
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
            eta = p_jh[0].m;
            for (size_t i=1;i<N;i++){
                reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r);
            }
            break;
    };
}

void reb_integrator_whfast_com_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    p_jh[0].x += _dt*p_jh[0].vx;
    p_jh[0].y += _dt*p_jh[0].vy;
    p_jh[0].z += _dt*p_jh[0].vz;
    // Only WHFast supports variational equations
    for (size_t v=0;v<r->N_var_config;v++){
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        struct reb_variational_configuration const vc = r->var_config[v];
        struct reb_particle* p_jh_var = whfast->p_jh_var;
        p_jh_var[vc.index].x += _dt*p_jh_var[vc.index].vx;
        p_jh_var[vc.index].y += _dt*p_jh_var[vc.index].vy;
        p_jh_var[vc.index].z += _dt*p_jh_var[vc.index].vz;
    }
}

static void reb_whfast_corrector_Z(struct reb_simulation* r, const double a, const double b){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
    switch (whfast->coordinates){
        case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);
            reb_transformations_jacobi_to_inertial_pos(particles, whfast->p_jh, particles, N, N_active);
            for (size_t v=0;v<r->N_var_config;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_jacobi_to_inertial_pos(r->particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
            }
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, -b);
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -2.*a);
            reb_transformations_jacobi_to_inertial_pos(particles, whfast->p_jh, particles, N, N_active);
            for (size_t v=0;v<r->N_var_config;v++){
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_jacobi_to_inertial_pos(r->particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
            }
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, b);
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);
            reb_transformations_barycentric_to_inertial_pos(particles, whfast->p_jh, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, -b);
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -2.*a);
            reb_transformations_barycentric_to_inertial_pos(particles, whfast->p_jh, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, b);
            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
        case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            reb_simulation_error(r, "Coordinate system not supported.");
            break;
    }
}

static void reb_whfast_apply_corrector(struct reb_simulation* r, double inv, int order){
    const double dt = r->dt;
    if (order==3){
        // Third order corrector
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_31*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_31*dt);
    }
    if (order==5){
        // Fifth order corrector
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_51*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_52*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_52*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_51*dt);
    }
    if (order==7){
        // Seventh order corrector
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_71*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_72*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_73*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_73*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_72*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_71*dt);
    }
    if (order==11){
        // Eleventh order corrector
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_5*dt,-inv*reb_whfast_corrector_b_111*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_4*dt,-inv*reb_whfast_corrector_b_112*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_113*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_114*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_115*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_115*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_114*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_113*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_4*dt,inv*reb_whfast_corrector_b_112*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_5*dt,inv*reb_whfast_corrector_b_111*dt);
    }
    if (order==17){
        // Seventeenth order corrector
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_8*dt,-inv*reb_whfast_corrector_b_171*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_7*dt,-inv*reb_whfast_corrector_b_172*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_6*dt,-inv*reb_whfast_corrector_b_173*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_5*dt,-inv*reb_whfast_corrector_b_174*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_4*dt,-inv*reb_whfast_corrector_b_175*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_176*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_177*dt);
        reb_whfast_corrector_Z(r, -reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_178*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_178*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_177*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_176*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_4*dt,inv*reb_whfast_corrector_b_175*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_5*dt,inv*reb_whfast_corrector_b_174*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_6*dt,inv*reb_whfast_corrector_b_173*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_7*dt,inv*reb_whfast_corrector_b_172*dt);
        reb_whfast_corrector_Z(r, reb_whfast_corrector_a_8*dt,inv*reb_whfast_corrector_b_171*dt);
    }
}

static void reb_whfast_operator_C(struct reb_simulation* const r, double a, double b){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);   

    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
    reb_transformations_jacobi_to_inertial_pos(particles, whfast->p_jh, particles, N, N_active); 
    // Note: variational particles not implemented.
    reb_simulation_update_acceleration(r);
    reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, b);

    reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -a);   
}

static void reb_whfast_operator_Y(struct reb_simulation* const r, double a, double b){
    reb_whfast_operator_C(r, a, b); 
    reb_whfast_operator_C(r, -a, -b); 
}
static void reb_whfast_operator_U(struct reb_simulation* const r, double a, double b){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, a);   
    reb_whfast_operator_Y(r, a, b); 
    reb_whfast_operator_Y(r, a, -b); 
    reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -a);   
}
static void reb_whfast_apply_corrector2(struct reb_simulation* r, double inv){
    double a = 0.5 * inv * r->dt;
    double b = reb_whfast_corrector2_b * inv * r->dt;
    reb_whfast_operator_U(r, a, b); 
    reb_whfast_operator_U(r, -a, b);
}

void reb_integrator_whfast_calculate_jerk(struct reb_simulation* r, struct reb_particle* jerk){
    // Assume particles.a calculated.
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    double Rjx = 0.; // com
    double Rjy = 0.;
    double Rjz = 0.;
    double Mj = 0.;
    double Ajx = 0.; // sort of Jacobi acceleration
    double Ajy = 0.;
    double Ajz = 0.;
    for (size_t j=0; j<N; j++){
        jerk[j].ax = 0; 
        jerk[j].ay = 0; 
        jerk[j].az = 0; 
        for (size_t i=0; i<j+1; i++){
            //////////////////
            // Jacobi Term
            // Note: ignoring j==1 term here and below as they cancel
            if (j>1){
                double dQkrj = Mj;
                if (i<j){
                    dQkrj = -particles[j].m;
                }
                const double Qkx = particles[j].x - Rjx/Mj;
                const double Qky = particles[j].y - Rjy/Mj;
                const double Qkz = particles[j].z - Rjz/Mj;
                const double dax = particles[j].ax - Ajx/Mj;
                const double day = particles[j].ay - Ajy/Mj;
                const double daz = particles[j].az - Ajz/Mj;

                const double dr = sqrt(Qkx*Qkx + Qky*Qky + Qkz*Qkz);

                const double prefact2 = G*dQkrj /(dr*dr*dr);
                jerk[i].ax    += prefact2*dax;
                jerk[i].ay    += prefact2*day;
                jerk[i].az    += prefact2*daz;

                const double alphasum = dax*Qkx + day*Qky + daz*Qkz;
                const double prefact1 = 3.*alphasum*prefact2/(dr*dr);
                jerk[i].ax    -= prefact1*Qkx; 
                jerk[i].ay    -= prefact1*Qky;
                jerk[i].az    -= prefact1*Qkz; 
            }
            /////////////////
            // Direct Term
            // Note: ignoring i==0 && j==1 term here and above as they cancel
            if (j!=i && (i!=0 || j!=1)){
                const double dx = particles[j].x - particles[i].x; 
                const double dy = particles[j].y - particles[i].y; 
                const double dz = particles[j].z - particles[i].z; 

                const double dax = particles[j].ax - particles[i].ax; 
                const double day = particles[j].ay - particles[i].ay; 
                const double daz = particles[j].az - particles[i].az; 

                const double dr = sqrt(dx*dx + dy*dy + dz*dz);
                const double alphasum = dax*dx+day*dy+daz*dz;
                const double prefact2 = G /(dr*dr*dr);
                const double prefact2i = prefact2*particles[i].m;
                const double prefact2j = prefact2*particles[j].m;
                jerk[j].ax    -= dax*prefact2i;
                jerk[j].ay    -= day*prefact2i;
                jerk[j].az    -= daz*prefact2i;
                jerk[i].ax    += dax*prefact2j;
                jerk[i].ay    += day*prefact2j;
                jerk[i].az    += daz*prefact2j;
                const double prefact1 = 3.*alphasum*prefact2 /(dr*dr);
                const double prefact1i = prefact1*particles[i].m;
                const double prefact1j = prefact1*particles[j].m;
                jerk[j].ax    += dx*prefact1i;
                jerk[j].ay    += dy*prefact1i;
                jerk[j].az    += dz*prefact1i;
                jerk[i].ax    -= dx*prefact1j;
                jerk[i].ay    -= dy*prefact1j;
                jerk[i].az    -= dz*prefact1j;
            }
        }
        Ajx += particles[j].ax*particles[j].m;
        Ajy += particles[j].ay*particles[j].m;
        Ajz += particles[j].az*particles[j].m;
        Rjx += particles[j].x*particles[j].m;
        Rjy += particles[j].y*particles[j].m;
        Rjz += particles[j].z*particles[j].m;
        Mj += particles[j].m;
    }
}


int reb_integrator_whfast_init(struct reb_simulation* const r, struct reb_integrator_whfast_state* whfast){
    if (r->N_var){
        if (whfast->coordinates!=REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI){
            reb_simulation_error(r, "Variational particles are only compatible with Jacobi coordinates.");
            return 1; // Error
        }
        if (whfast->kernel != REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT){
            reb_simulation_error(r, "Variational particles are only compatible with the standard kernel.");
            return 1; // Error
        }
        if (whfast->corrector2){
            reb_simulation_error(r, "Variational particles not compatible with 2nd corrector.");
            return 1; // Error
        }
        if (whfast->safe_mode==0&&r->calculate_megno){
            reb_simulation_error(r, "MEGNO is not compatible with WHFast's safe_mode=0.");
            return 1; // Error
        }
        if (whfast->N_allocated_var != r->N_var){
            whfast->N_allocated_var = r->N_var;
            whfast->p_jh_var = realloc(whfast->p_jh_var, sizeof(struct reb_particle)*r->N_var);
        }
        for (size_t v=0;v<r->N_var_config;v++){
            struct reb_variational_configuration const vc = r->var_config[v];
            if (vc.order!=1){
                reb_simulation_error(r, "WHFast only supports first order variational equations.");
                return 1; // Error
            }
            if (vc.testparticle>=0){
                reb_simulation_error(r, "Test particle variations not supported with WHFast. Use IAS15.");
                return 1; // Error
            }
        }
    }
    if (whfast->kernel!= REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT && whfast->coordinates!=REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI){
        reb_simulation_error(r, "Non-standard kernel requires Jacobi coordinates.");
        return 1; // Error
    }
    if (whfast->kernel>3){
        reb_simulation_error(r, "Kernel method must be 0 (default), 1 (exact modified kick), 2 (composition kernel), or 3 (lazy implementer's modified kick). ");
        return 1; // Error
    }
    if (whfast->corrector!=0 && (whfast->coordinates!=REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI && whfast->coordinates!=REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC) ){
        reb_simulation_error(r, "Symplectic correctors are only compatible with Jacobi and Barycentric coordinates.");
        return 1; // Error
    }
    if (whfast->corrector!=0 && whfast->corrector!=3 && whfast->corrector!=5  && whfast->corrector!=7 && whfast->corrector!=11 && whfast->corrector!=17 ){
        reb_simulation_error(r, "First symplectic correctors are only available in the following orders: 0, 3, 5, 7, 11, 17.");
        return 1; // Error
    }
    if (whfast->keep_unsynchronized==1 && whfast->safe_mode==1){
        reb_simulation_error(r, "whfast->keep_unsynchronized == 1 is not compatible with safe_mode. Must set whfast->safe_mode = 0.");
    }
    if (whfast->kernel == REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK || whfast->kernel == REB_INTEGRATOR_WHFAST_KERNEL_LAZY){ 
        r->gravity = REB_GRAVITY_JACOBI;
    }else{
        if (whfast->coordinates==REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI){
            r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1;
        }else if (whfast->coordinates==REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC){
            r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
        }else{
            r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_INVOLVING_0;
        }
    }
    const size_t N = r->N;
    if (whfast->N_allocated != N){
        whfast->N_allocated = N;
        whfast->p_jh = realloc(whfast->p_jh,sizeof(struct reb_particle)*N);
        whfast->recalculate_coordinates_this_timestep = 1;
    }
    return 0;
}

void reb_integrator_whfast_from_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;

    switch (coordinates){
        case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
            reb_transformations_inertial_to_jacobi_posvel(particles, p_jh, particles, N, N_active);
            // Only WHFast supports variational equations
            for (size_t v=0;v<r->N_var_config;v++){
                struct reb_integrator_whfast_state* whfast = r->integrator.state;
                struct reb_variational_configuration const vc = r->var_config[v];
                reb_transformations_inertial_to_jacobi_posvel(r->particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
            }
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
            reb_transformations_inertial_to_democraticheliocentric_posvel(particles, p_jh, N, N_active);
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
            reb_transformations_inertial_to_whds_posvel(particles, p_jh, N, N_active);
            break;
        case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
            reb_transformations_inertial_to_barycentric_posvel(particles, p_jh, N, N_active);
            break;
    };
}

void reb_integrator_whfast_to_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;

    // Prepare coordinates for KICK step
    if (r->force_is_velocity_dependent){
        switch (coordinates){
            case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(particles, p_jh, particles, N, N_active);
                for (size_t v=0;v<r->N_var_config;v++){
                    // Only WHFast supports variational equations
                    struct reb_integrator_whfast_state* whfast = r->integrator.state;
                    struct reb_variational_configuration const vc = r->var_config[v];
                    reb_transformations_jacobi_to_inertial_posvel(r->particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
                }
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
                reb_transformations_barycentric_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
        };
    }else{
        switch (coordinates){
            case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(particles, p_jh, particles, N, N_active);
                for (size_t v=0;v<r->N_var_config;v++){
                    struct reb_integrator_whfast_state* whfast = r->integrator.state;
                    struct reb_variational_configuration const vc = r->var_config[v];
                    reb_transformations_jacobi_to_inertial_pos(r->particles_var+vc.index, whfast->p_jh_var+vc.index, particles, N, N_active);
                }
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
                reb_transformations_barycentric_to_inertial_posvel(particles, p_jh, N, N_active);
                break;
        };
    }
}

void reb_integrator_whfast_debug_operator_kepler(struct reb_simulation* const r,double dt){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    if (reb_integrator_whfast_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }
    reb_integrator_whfast_from_inertial(r, whfast->p_jh, whfast->coordinates);
    reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, dt);    // half timestep
    reb_integrator_whfast_com_step(r, whfast->p_jh, dt);
    reb_integrator_whfast_to_inertial(r, whfast->p_jh, whfast->coordinates);
}

void reb_integrator_whfast_debug_operator_interaction(struct reb_simulation* const r,double dt){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    if (reb_integrator_whfast_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }
    reb_integrator_whfast_from_inertial(r, whfast->p_jh, whfast->coordinates);
    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1;
    reb_simulation_update_acceleration(r);
    reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);
    reb_integrator_whfast_to_inertial(r, whfast->p_jh, whfast->coordinates);
}

void reb_integrator_whfast_synchronize(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_state* whfast = state;
    if (reb_integrator_whfast_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }
    if (r->is_synchronized == 0){
        const size_t N = r->N;
        const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
        struct reb_particle* sync_pj  = NULL;
        if (whfast->keep_unsynchronized){
            sync_pj = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(sync_pj,whfast->p_jh,r->N*sizeof(struct reb_particle));
        }
        switch (whfast->kernel){
            case REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT: 
            case REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK: 
            case REB_INTEGRATOR_WHFAST_KERNEL_LAZY: 
                reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, r->dt/2.);    
                reb_integrator_whfast_com_step(r, whfast->p_jh, r->dt/2.);
                break;
            case REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION:
                reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, 3.*r->dt/8.);   
                reb_integrator_whfast_com_step(r, whfast->p_jh, 3.*r->dt/8.);
                break;
            default:
                reb_simulation_error(r, "WHFast kernel not implemented.");
                return;
        };
        if (whfast->corrector2){
            reb_whfast_apply_corrector2(r, -1.);
        }
        if (whfast->corrector){
            reb_whfast_apply_corrector(r, -1., whfast->corrector);
        }
        switch (whfast->coordinates){
            case REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI:
                reb_transformations_jacobi_to_inertial_posvel(r->particles, whfast->p_jh, r->particles, N, N_active);
                for (size_t v=0;v<r->N_var_config;v++){
                    struct reb_variational_configuration const vc = r->var_config[v];
                    reb_transformations_jacobi_to_inertial_posvel(r->particles_var+vc.index, whfast->p_jh_var+vc.index, r->particles, N, N_active);
                }
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC:
                reb_transformations_democraticheliocentric_to_inertial_posvel(r->particles, whfast->p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_WHDS:
                reb_transformations_whds_to_inertial_posvel(r->particles, whfast->p_jh, N, N_active);
                break;
            case REB_INTEGRATOR_WHFAST_COORDINATES_BARYCENTRIC:
                reb_transformations_barycentric_to_inertial_posvel(r->particles, whfast->p_jh, N, N_active);
                break;
        };
        if (whfast->keep_unsynchronized){
            memcpy(whfast->p_jh,sync_pj,r->N*sizeof(struct reb_particle));
            free(sync_pj);
        }else{
            r->is_synchronized = 1;
        }
    }
}

void reb_integrator_whfast_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_state* whfast = state;
    struct reb_particle* restrict const particles = r->particles;
    const double dt = r->dt;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?N:r->N_active;
    if (reb_integrator_whfast_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }
    struct reb_particle* const p_jh = whfast->p_jh;

    // Only recalculate Jacobi coordinates if needed
    if (whfast->safe_mode || whfast->recalculate_coordinates_this_timestep){
        if (r->is_synchronized==0){
            reb_integrator_whfast_synchronize(r, whfast);
            if (whfast->recalculate_coordinates_but_not_synchronized_warning==0){
                reb_simulation_warning(r,"Recalculating coordinates but pos/vel were not synchronized before.");
                whfast->recalculate_coordinates_but_not_synchronized_warning++;
            }
        }
        reb_integrator_whfast_from_inertial(r, whfast->p_jh, whfast->coordinates);
        whfast->recalculate_coordinates_this_timestep = 0;
    }
    if (r->is_synchronized){
        // First half DRIFT step
        if (whfast->corrector){
            reb_whfast_apply_corrector(r, 1., whfast->corrector);
        }
        if (whfast->corrector2){
            reb_whfast_apply_corrector2(r, 1.);
        }
        switch (whfast->kernel){
            case REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT: 
            case REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK: 
            case REB_INTEGRATOR_WHFAST_KERNEL_LAZY: 
                reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, r->dt/2.);    
                reb_integrator_whfast_com_step(r, whfast->p_jh, r->dt/2.);
                break;
            case REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION:
                reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, 5.*r->dt/8.);   
                reb_integrator_whfast_com_step(r, whfast->p_jh, 5.*r->dt/8.);
                break;
            default:
                reb_simulation_error(r, "WHFast kernel not implemented.");
                return;
        };
    }else{
        // Combined DRIFT step
        reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, r->dt);    // full timestep
        reb_integrator_whfast_com_step(r, whfast->p_jh, r->dt);
    }
    reb_integrator_whfast_jump_step(r, whfast, r->dt/2.);

    reb_integrator_whfast_to_inertial(r, whfast->p_jh, whfast->coordinates);

    r->t+=dt/2.;

    reb_simulation_update_acceleration(r);

    switch (whfast->kernel){
        case REB_INTEGRATOR_WHFAST_KERNEL_DEFAULT: 
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);

            reb_integrator_whfast_jump_step(r, whfast, dt/2.);
            break;
        case REB_INTEGRATOR_WHFAST_KERNEL_MODIFIEDKICK: 
            // p_jh used as a temporary buffer for "jerk"
            reb_integrator_whfast_calculate_jerk(r, whfast->p_jh);
            for (size_t i=0; i<N; i++){
                const double prefact = dt*dt/12.;
                particles[i].ax += prefact*p_jh[i].ax; 
                particles[i].ay += prefact*p_jh[i].ay; 
                particles[i].az += prefact*p_jh[i].az; 
            }
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);
            break;
        case REB_INTEGRATOR_WHFAST_KERNEL_COMPOSITION:
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, -dt/6.);

            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -dt/4.);   
            reb_integrator_whfast_com_step(r, whfast->p_jh, -dt/4.);

            reb_transformations_jacobi_to_inertial_pos(particles, p_jh, particles, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt/6.);

            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, dt/8.);   
            reb_integrator_whfast_com_step(r, whfast->p_jh, dt/8.);

            reb_transformations_jacobi_to_inertial_pos(particles, p_jh, particles, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);

            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, -dt/8.);   
            reb_integrator_whfast_com_step(r, whfast->p_jh, -dt/8.);

            reb_transformations_jacobi_to_inertial_pos(particles, p_jh, particles, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, -dt/6.);

            reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, dt/4.);   
            reb_integrator_whfast_com_step(r, whfast->p_jh, dt/4.);

            reb_transformations_jacobi_to_inertial_pos(particles, p_jh, particles, N, N_active);
            reb_simulation_update_acceleration(r);
            reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt/6.);
            break;
        case REB_INTEGRATOR_WHFAST_KERNEL_LAZY: 
            {
                // Need temporary array to store old positions
                if (whfast->N_allocated_temp != N){
                    whfast->N_allocated_temp = N;
                    whfast->p_temp = realloc(whfast->p_temp,sizeof(struct reb_particle)*N);
                }
                struct reb_particle* p_temp = whfast->p_temp;

                // Calculate normal kick
                // Accelerations already calculated
                reb_transformations_inertial_to_jacobi_acc(r->particles, p_jh, r->particles, N, N_active);

                // make copy of original positions
                memcpy(p_temp,p_jh,r->N*sizeof(struct reb_particle));

                // WHT Eq 10.6
                for (size_t i=1;i<N;i++){
                    const double prefac1 = dt*dt/12.; 
                    p_jh[i].x += prefac1 * p_temp[i].ax;
                    p_jh[i].y += prefac1 * p_temp[i].ay;
                    p_jh[i].z += prefac1 * p_temp[i].az;
                }

                // recalculate kick 
                reb_transformations_jacobi_to_inertial_pos(particles, p_jh, particles, N, N_active);
                reb_simulation_update_acceleration(r);
                reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);

                for (size_t i=1;i<N;i++){
                    // reset positions
                    p_jh[i].x = p_temp[i].x;
                    p_jh[i].y = p_temp[i].y;
                    p_jh[i].z = p_temp[i].z;
                }
            }
            break;
        default:
            return;
    };

    r->is_synchronized = 0;
    if (whfast->safe_mode){
        reb_integrator_whfast_synchronize(r, whfast);
    }

    r->t+=r->dt/2.;
    r->dt_last_done = r->dt;

    if (r->calculate_megno){
        // Need to have x,v,a synchronized to calculate ddot/d for MEGNO. 
        r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE; // Need all terms.
        reb_gravity_basic_calculate_acceleration_var(r);
        r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_INVOLVING_0;

        double dY = r->dt * 2. * (r->t-r->megno_initial_t) * reb_tools_megno_deltad_delta(r);
        reb_tools_megno_update(r, dY, dt);
    }
}

