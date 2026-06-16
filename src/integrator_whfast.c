/**
 * integrator_whfast.c: The Wisdom-Holman integrator WHFast
 * 
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
#include "transformations.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast.h"
#include "binarydata.h"
#include <gmp.h>
#include <mpfr.h>
#include <fenv.h>

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

void* reb_integrator_whfast_create();
void reb_integrator_whfast_free(void* state);
void reb_integrator_whfast_synchronize(struct reb_simulation* const r, void* state);
void reb_integrator_whfast_step(struct reb_simulation* const r, void* state);
const struct reb_binarydata_field_descriptor reb_integrator_whfast_field_descriptor_list[];

const struct reb_integrator reb_integrator_whfast = {
    .documentation =
    "WHFast is an implementation of the symplectic [Wisdom & Holman (1991)] integrator. " 
    "It is the best choice for systems in which there is a dominant central object and perturbations to the Keplerian orbits are small. "
    "It supports first and second symplectic correctors as well as the kernel method of [Wisdom et al. (1996)] with various different kernels. "
    "The basic implementation of WHFast is described in detail in [Rein & Tamayo (2015)]. "
    "The higher order aspects of it are described in [Rein, Tamayo & Brown (2019)]. "
    "WHFast also supports first order variational equations which can be used in chaos estimators. See [Rein & Tamayo (2016)]. "
    "The user can choose between Jacobi, Democratic Heliocentric, WHDS, and barycentric coordinates. "
    "Because WHFast is not an adaptive integrator, the user needs to set an appropriaye timestep. "
    "Typically, this should be a small fraction (a few percent) of the smallest dynamical timescale in the problem. "
    "\n\n"
    "[Wisdom & Holman (1991)]: https://ui.adsabs.harvard.edu/abs/1991AJ....102.1528W/abstract\n"
    "[Wisdom et al. (1996)]: https://ui.adsabs.harvard.edu/abs/1996FIC....10..217W/abstract\n"
    "[Rein & Tamayo (2015)]: https://ui.adsabs.harvard.edu/abs/2015MNRAS.452..376R/abstract\n"
    "[Rein, Tamayo & Brown (2019)]: https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R/abstract\n" 
    "[Rein & Tamayo (2016)]: https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2275R/abstract\n" 
    "[Rein & Tamayo (2016)]: https://www.youtube.com/watch?v=ttLUhtNj1Lc\n" 
    ,
    .step = reb_integrator_whfast_step,
    .synchronize = reb_integrator_whfast_synchronize,
    .create = reb_integrator_whfast_create,
    .free = reb_integrator_whfast_free,
    .field_descriptor_list = reb_integrator_whfast_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast_field_descriptor_list[] = {
    { "This flag sets the order of the symplectic corrector in the WHFast integrator. "
        "By default, the symplectic correctors are turned off (=0). For high "
        "accuracy simulation set this value to 11 or 17. For more details read "
        "Rein and Tamayo (2015). ",
        REB_UINT,       "corrector",          offsetof(struct reb_integrator_whfast_state, corrector), 0, 0, 0},
    { "If this flag is set to 1 (default) particle positions and velocities are "
        "always synchronized and particles can be modified between timesteps. "
        "If this flag is set to 0, the speed and accuracy of WHFast improve. "
        "However, one needs to make sure to call synchronize before an output is "
        "required or before particles are modified. Read the iPython tutorial "
        "on advanced WHFast usage to learn more.",
        REB_UINT,        "safe_mode",          offsetof(struct reb_integrator_whfast_state, safe_mode), 0, 0, 0},
    { "This option chooses the internal coordinate system that WHFast is using for "
        "splitting the Keplerian from the Interaction and Jump parts. "
        "By default, it uses JACOBI coordinates. See Rein & Tamayo 2019 and "
        "Hernandez & Dehnen (2017) for more information.",
        REB_INT,         "coordinates",        offsetof(struct reb_integrator_whfast_state, coordinates), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_WHFAST_COORDINATES) },
    { "This flag chooses the second correctors (C2 of Wisdom et al 1996)."
        " By default, the second symplectic correctors are turned off (=0). "
        " Set to 1 to turn them on. ",
        REB_UINT,        "corrector2",         offsetof(struct reb_integrator_whfast_state, corrector2), 0, 0, 0},
    { "This flag chooses the kernel. The default is DEFAULT which corresponds "
        "to the normal 2nd order WH kernel (i.e. a standard kick step). "
        "See Rein, Tamayo & Brown 2019 for details and references. ",
        REB_INT,         "kernel",             offsetof(struct reb_integrator_whfast_state, kernel), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_WHFAST_KERNEL)},
    { "By default the keep_unsynchronized flag is 0. If set to 1 synchronization of the "
        "simulation is done on a copy of the particle data. This allows "
        "the simulation to continue integrating as if the simulation "
        "were never synchronized. This allows for bit-wise reproducibility "
        "in long term simulations.",
        REB_UINT,        "keep_unsynchronized",offsetof(struct reb_integrator_whfast_state, keep_unsynchronized), 0, 0, 0},
    // Internal variables
    { "", REB_POINTER,     "p_jh",               offsetof(struct reb_integrator_whfast_state, p_jh), offsetof(struct reb_integrator_whfast_state, N_allocated), sizeof(struct reb_particle), 0},
    { "", REB_POINTER,     "p_jh_var",           offsetof(struct reb_integrator_whfast_state, p_jh_var), offsetof(struct reb_integrator_whfast_state, N_allocated_var), sizeof(struct reb_particle), 0},
    { 0 }, // Null terminated list
};


void* reb_integrator_whfast_create(){
    // Allocate memory and set default parameters.
    struct reb_integrator_whfast_state* whfast = calloc(sizeof(struct reb_integrator_whfast_state),1);
    whfast->coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_JACOBI;
    whfast->safe_mode = 1;
    return whfast;
}

void reb_integrator_whfast_free(void* p){
    struct reb_integrator_whfast_state* whfast = p;
    free(whfast->p_jh);
    free(whfast->p_jh_var);
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

static inline double fastabs(double x){
    return (x > 0.) ? x : -x;
}


void mpfr_invfactorial(mpfr_t fac, unsigned int i){
    mpfr_set_d(fac, 1.0, MPFR_RNDN);
    for(size_t l=2;l<=i;l++){
        mpfr_mul_ui(fac, fac, l, MPFR_RNDN);
    }
    mpfr_ui_div(fac, 1, fac, MPFR_RNDN);
}

mpfr_t ifac, c_odd, c_even, z, tmp, cs0, cs1, cs2, cs3, mX, mri, mf, mg, mfd, mgd, mr0i, mr0, mbeta, mzeta0, mv2, meta0, mx, my, mz, mvx, mvy, mvz, tmpx, tmpy, tmpz;
void init_mpfr(){
    mpfr_init2(mX,300);
    mpfr_init2(mr0i,300);
    mpfr_init2(z,300);
    mpfr_init2(mx,300);
    mpfr_init2(my,300);
    mpfr_init2(mz,300);
    mpfr_init2(tmpx,300);
    mpfr_init2(tmpy,300);
    mpfr_init2(tmpz,300);
    mpfr_init2(mvx,300);
    mpfr_init2(mvy,300);
    mpfr_init2(mvz,300);
    mpfr_init2(mr0,300);
    mpfr_init2(mbeta,300);
    mpfr_init2(meta0,300);
    mpfr_init2(mzeta0,300);
    mpfr_init2(mv2,300);
    mpfr_init2(mf,300);
    mpfr_init2(mg,300);
    mpfr_init2(mgd,300);
    mpfr_init2(mfd,300);
    mpfr_init2(mri,300);
    mpfr_init2(cs0,300);
    mpfr_init2(cs1,300);
    mpfr_init2(cs2,300);
    mpfr_init2(cs3,300);
    mpfr_init2(tmp,300);
    mpfr_init2(ifac,300);
    mpfr_init2(c_odd,300);
    mpfr_init2(c_even,300);
}
static void stumpff_cs3(double *restrict cs, double _z) {
    (void)cs;
    (void)_z;
    unsigned int n = 0;
    while(1){
        mpfr_abs(tmp, z, MPFR_RNDN);
        if (mpfr_cmp_d(tmp, 0.1)<0.0) break;
        mpfr_div_ui(z, z, 4, MPFR_RNDN);
        n++;
    }
    const int nmax = 15;
    mpfr_invfactorial(c_odd, nmax);
    mpfr_mul_ui(c_even, c_odd, nmax, MPFR_RNDN);
    mpfr_set(ifac, c_even, MPFR_RNDN);

    for(int np=nmax-2;np>=3;np-=2){
        mpfr_mul_ui(ifac, ifac, np+1, MPFR_RNDN);
        mpfr_mul(c_odd, c_odd, z, MPFR_RNDN);
        mpfr_sub(c_odd, ifac, c_odd, MPFR_RNDN);
        mpfr_mul_ui(ifac, ifac, np, MPFR_RNDN);
        mpfr_mul(c_even, c_even, z, MPFR_RNDN);
        mpfr_sub(c_even, ifac, c_even, MPFR_RNDN);
    }
    mpfr_set(cs3,c_odd, MPFR_RNDN);
    mpfr_set(cs2,c_even, MPFR_RNDN);
    mpfr_mul(c_odd, c_odd, z, MPFR_RNDN);
    mpfr_mul(c_even, c_even, z, MPFR_RNDN);
    mpfr_ui_sub(c_odd, 1, c_odd, MPFR_RNDN);
    mpfr_ui_sub(c_even, 1, c_even, MPFR_RNDN);
    mpfr_set(cs1,c_odd, MPFR_RNDN);
    mpfr_set(cs0,c_even, MPFR_RNDN);
    for (;n>0;n--){ 
        mpfr_mul(cs3, cs3, cs0, MPFR_RNDN);
        mpfr_add(cs3, cs3, cs2, MPFR_RNDN);
        mpfr_div_ui(cs3, cs3, 4, MPFR_RNDN);
        mpfr_mul(cs2, cs1, cs1, MPFR_RNDN);
        mpfr_div_ui(cs2, cs2, 2, MPFR_RNDN);
        mpfr_mul(cs1, cs1, cs0, MPFR_RNDN);
        mpfr_mul(cs0, cs1, cs0, MPFR_RNDN);
        mpfr_mul_ui(cs0, cs0, 2, MPFR_RNDN);
        mpfr_sub_ui(cs0, cs0, 1, MPFR_RNDN);
    }
}

//static void stiefel_Gs(double *restrict Gs, double beta, double X) {
//    // Only variations use this
//    double X2 = X*X;
//    stumpff_cs(Gs, beta*X2);
//    Gs[1] *= X; 
//    Gs[2] *= X2; 
//    double _pow = X2*X;
//    Gs[3] *= _pow; 
//    _pow *= X;
//    Gs[4] *= _pow; 
//    _pow *= X;
//    Gs[5] *= _pow; 
//    return;
//}

static void stiefel_Gs3() {
    //mpfr_set_d(mX, X, MPFR_RNDN);
    mpfr_mul(z, mX, mX, MPFR_RNDN);
    mpfr_mul(z, z, mbeta, MPFR_RNDN);
    stumpff_cs3(NULL, 0);
    mpfr_mul(cs1, cs1, mX, MPFR_RNDN);
    mpfr_mul(cs2, cs2, mX, MPFR_RNDN);
    mpfr_mul(cs2, cs2, mX, MPFR_RNDN);
    mpfr_mul(cs3, cs3, mX, MPFR_RNDN);
    mpfr_mul(cs3, cs3, mX, MPFR_RNDN);
    mpfr_mul(cs3, cs3, mX, MPFR_RNDN);
    return;
}

#define WHFAST_NMAX_QUART 64    ///< Maximum number of iterations for quartic solver
#define WHFAST_NMAX_NEWT  10    ///< Maximum number of iterations for Newton's method
                                // Keplerian motion for one planet                       
                                // r only needed for variational particles and warning. Can be NULL.
void reb_integrator_whfast_kepler_solver(struct reb_particle* const restrict p, double mu, double dt, const struct reb_simulation* const r){
    const struct reb_particle p1 = *p; // Copy of particle

    if (r->t==0){
        mpfr_set_d(mx, p1.x, MPFR_RNDN);
        mpfr_set_d(my, p1.y, MPFR_RNDN);
        mpfr_set_d(mz, p1.z, MPFR_RNDN);
        mpfr_set_d(mvx, p1.vx, MPFR_RNDN);
        mpfr_set_d(mvy, p1.vy, MPFR_RNDN);
        mpfr_set_d(mvz, p1.vz, MPFR_RNDN);
    }

    mpfr_set(mr0, mx, MPFR_RNDN);
    mpfr_mul(mr0, mr0, mr0, MPFR_RNDN);
    mpfr_set(tmp, my, MPFR_RNDN);
    mpfr_fma(mr0, tmp, tmp, mr0, MPFR_RNDN);
    mpfr_set(tmp, mz, MPFR_RNDN);
    mpfr_fma(mr0, tmp, tmp, mr0, MPFR_RNDN);
    mpfr_sqrt(mr0, mr0, MPFR_RNDN);
    mpfr_ui_div(mr0i, 1, mr0, MPFR_RNDN);

    const double r0i = mpfr_get_d(mr0i, MPFR_RNDN);
    
    mpfr_set(mv2, mvx, MPFR_RNDN);
    mpfr_mul(mv2, mv2, mv2, MPFR_RNDN);
    mpfr_set(tmp, mvy, MPFR_RNDN);
    mpfr_fma(mv2, tmp, tmp, mv2, MPFR_RNDN);
    mpfr_set(tmp, mvz, MPFR_RNDN);
    mpfr_fma(mv2, tmp, tmp, mv2, MPFR_RNDN);
    
    mpfr_mul_ui(mbeta, mr0i, 2, MPFR_RNDN);
    mpfr_sub(mbeta, mbeta, mv2, MPFR_RNDN);

    mpfr_set(tmp, mx, MPFR_RNDN);
    mpfr_mul(meta0, tmp, mvx, MPFR_RNDN);
    mpfr_set(tmp, my, MPFR_RNDN);
    mpfr_mul(tmp, tmp, mvy, MPFR_RNDN);
    mpfr_add(meta0, tmp,meta0, MPFR_RNDN);
    mpfr_set(tmp, mz, MPFR_RNDN);
    mpfr_mul(tmp, tmp, mvz, MPFR_RNDN);
    mpfr_add(meta0, tmp,meta0, MPFR_RNDN);

    mpfr_mul(tmp, mbeta, mr0, MPFR_RNDN);
    mpfr_ui_sub(mzeta0, 1.0, tmp, MPFR_RNDN);



    const double beta = mpfr_get_d(mbeta,MPFR_RNDN);
    const double eta0 = mpfr_get_d(meta0, MPFR_RNDN);
    double X;
    double invperiod=0;  // only used for beta>0. Set to 0 only to suppress compiler warnings.

    if (beta>0.){
        // Elliptic orbit
        const double sqrt_beta = sqrt(beta);
        invperiod = sqrt_beta*beta/(2.*M_PI*mu);
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
        printf("hyperbolic\n");
        exit(1);
    }
    mpfr_set_d(mX, X, MPFR_RNDN);

    for (int n_hg=1;n_hg<WHFAST_NMAX_NEWT;n_hg++){
        stiefel_Gs3();
        mpfr_mul(mri, cs1, meta0, MPFR_RNDN);
        mpfr_mul(tmp, cs2, mzeta0, MPFR_RNDN);
        mpfr_add(tmp, mri, tmp, MPFR_RNDN);
        mpfr_add(mri, tmp, mr0, MPFR_RNDN);
        mpfr_ui_div(mri, 1, mri, MPFR_RNDN);

        mpfr_mul(mX,tmp,mX, MPFR_RNDN);
        mpfr_mul(tmp,meta0,cs2, MPFR_RNDN);
        mpfr_sub(mX, mX,tmp, MPFR_RNDN);
        mpfr_mul(tmp,mzeta0,cs3, MPFR_RNDN);
        mpfr_sub(mX, mX,tmp, MPFR_RNDN);
        mpfr_add_d(mX, mX, dt, MPFR_RNDN);
        mpfr_mul(mX, mX, mri, MPFR_RNDN);
        //const double eta0Gs1zeta0Gs2 = eta0*Gs[1] + zeta0*Gs[2];
        //double ri = 1./(r0 + eta0Gs1zeta0Gs2);
        //X  = ri*(X*eta0Gs1zeta0Gs2-eta0*Gs[2]-zeta0*Gs[3]+dt);

        //if (X==oldX||X==oldX2){
        //    // Converged. Exit.
        //    converged = 1;
        //    break; 
        //}
    }

//    int original_mode = fegetround();
//    fesetround(FE_TOWARDZERO);
//    stiefel_Gs(Gs, beta, X);
    mpfr_mul(mri, cs1, meta0, MPFR_RNDN);
    mpfr_mul(tmp, cs2, mzeta0, MPFR_RNDN);
    mpfr_add(mri, mri, tmp, MPFR_RNDN);
    mpfr_add(mri, mri, mr0, MPFR_RNDN);
    mpfr_ui_div(mri, 1, mri, MPFR_RNDN);

    mpfr_mul(mf, cs2,mr0i, MPFR_RNDN);
    mpfr_neg(mf, mf, MPFR_RNDN);
    mpfr_d_sub(mg, dt,cs3, MPFR_RNDN);
    mpfr_mul(mfd, cs1,mri, MPFR_RNDN);
    mpfr_mul(mfd, cs1,mr0i, MPFR_RNDN);
    mpfr_neg(mfd, mfd, MPFR_RNDN);
    mpfr_mul(mfd, mfd,mri, MPFR_RNDN);
    mpfr_mul(mgd, cs2, mri, MPFR_RNDN);
    mpfr_neg(mgd, mgd, MPFR_RNDN);

    
    //f = fd*g-f*gd-gd;
    mpfr_fma(tmpx, mf, mx, mx, MPFR_RNDN);
    mpfr_fma(tmpy, mf, my, my, MPFR_RNDN);
    mpfr_fma(tmpz, mf, mz, mz, MPFR_RNDN);
    mpfr_fma(tmpx, mg, mvx, tmpx, MPFR_RNDN);
    mpfr_fma(tmpy, mg, mvy, tmpy, MPFR_RNDN);
    mpfr_fma(tmpz, mg, mvz, tmpz, MPFR_RNDN);
    
    mpfr_fma(mvx, mgd, mvx, mvx, MPFR_RNDN);
    mpfr_fma(mvy, mgd, mvy, mvy, MPFR_RNDN);
    mpfr_fma(mvz, mgd, mvz, mvz, MPFR_RNDN);
    mpfr_fma(mvx, mfd, mx, mvx, MPFR_RNDN);
    mpfr_fma(mvy, mfd, my, mvy, MPFR_RNDN);
    mpfr_fma(mvz, mfd, mz, mvz, MPFR_RNDN);

    mpfr_set(mx, tmpx, MPFR_RNDN);
    mpfr_set(my, tmpy, MPFR_RNDN);
    mpfr_set(mz, tmpz, MPFR_RNDN);

    p->x =  mpfr_get_d(mx, MPFR_RNDN);
    p->y =  mpfr_get_d(my, MPFR_RNDN);
    p->z =  mpfr_get_d(mz, MPFR_RNDN);
    p->vx = mpfr_get_d(mvx,MPFR_RNDN);
    p->vy = mpfr_get_d(mvy,MPFR_RNDN);
    p->vz = mpfr_get_d(mvz,MPFR_RNDN);
//    fesetround(original_mode);

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
        r->did_modify_particles = 1;
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

void reb_integrator_whfast_debug_operator_kepler(struct reb_simulation* const r,double dt, size_t Nsteps){
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    if (reb_integrator_whfast_init(r, whfast)){
        // Non recoverable error occurred.
        return;
    }
    reb_integrator_whfast_from_inertial(r, whfast->p_jh, whfast->coordinates);
    for (size_t i = 0;i<Nsteps;i++){
        reb_integrator_whfast_kepler_step(r, whfast->p_jh, whfast->coordinates, dt);
    }
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
    if (whfast->safe_mode || r->did_modify_particles){
        if (r->is_synchronized==0){
            reb_integrator_whfast_synchronize(r, whfast);
            if (whfast->recalculate_coordinates_but_not_synchronized_warning==0){
                reb_simulation_warning(r,"Particles were modified while simulation was not synchronized.");
                whfast->recalculate_coordinates_but_not_synchronized_warning++;
            }
        }
        reb_integrator_whfast_from_inertial(r, whfast->p_jh, whfast->coordinates);
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
                // Accelerations already calculated
                // WHT Eq 10.6
                for (size_t i=1;i<N;i++){
                    const double prefac1 = dt*dt/12.; 
                    r->particles[i].x += prefac1 * r->particles[i].ax;
                    r->particles[i].y += prefac1 * r->particles[i].ay;
                    r->particles[i].z += prefac1 * r->particles[i].az;
                }
                // Position will be overwritten in next jacobi_to_inertial transformation.

                // recalculate kick 
                reb_simulation_update_acceleration(r);
                reb_integrator_whfast_interaction_step(r, whfast->p_jh, whfast->coordinates, dt);

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

