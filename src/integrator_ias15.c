/**
 * @file    integrator_ias15.c
 * @brief   IAS15 integrator.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the IAS15 integration scheme.  
 * IAS stands for Integrator with Adaptive Step-size control, 15th 
 * order. This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * For more details see Rein & Spiegel 2014. Also see Everhart, 1985,
 * ASSL Vol. 115, IAU Colloq. 83, Dynamics of Comets, Their Origin 
 * and Evolution, 185 for the original implementation by Everhart.
 * Part of this code is based a function from the ORSE package.
 * See orsa.sourceforge.net for more details on their implementation.
 *
 * 
 * @section     LICENSE
 * Copyright (c) 2011-2012 Hanno Rein, Dave Spiegel.
 * Copyright (c) 2002-2004 Pasquale Tricarico.
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
// Uncomment the following line to generate numerical constants with extended precision.
//#define GENERATE_CONSTANTS
#ifdef GENERATE_CONSTANTS
#include <gmp.h>
void integrator_generate_constants(void);
#endif // GENERATE_CONSTANTS
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "tools.h"
#include "integrator.h"
#include "integrator_ias15.h"

/**
 * @brief Struct containing pointers to intermediate values
 */
struct reb_dpconst7 {
    double* const restrict p0;  ///< Temporary values at intermediate step 0 
    double* const restrict p1;  ///< Temporary values at intermediate step 1 
    double* const restrict p2;  ///< Temporary values at intermediate step 2 
    double* const restrict p3;  ///< Temporary values at intermediate step 3 
    double* const restrict p4;  ///< Temporary values at intermediate step 4 
    double* const restrict p5;  ///< Temporary values at intermediate step 5 
    double* const restrict p6;  ///< Temporary values at intermediate step 6 
};

// Helper functions for resetting the b and e coefficients
static void copybuffers(const struct reb_dpconst7 _a, const struct reb_dpconst7 _b, int N3);
static void predict_next_step(double ratio, int N3,  const struct reb_dpconst7 _e, const struct reb_dpconst7 _b, const struct reb_dpconst7 e, const struct reb_dpconst7 b);


/////////////////////////
//   Constants 

static const double safety_factor           = 0.25; /**< Maximum increase/deacrease of consecutve timesteps. */

// Gauss Radau spacings
static const double h[8]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626};
// Other constants
static const double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};
static const double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588};


// Weights for integration of a first order differential equation (Note: interval length = 2) 
static const double w[8] = {0.03125, 0.185358154802979278540728972807180754479812609, 0.304130620646785128975743291458180383736715043, 0.376517545389118556572129261157225608762708603, 0.391572167452493593082499533303669362149363727, 0.347014795634501068709955597003528601733139176, 0.249647901329864963257869294715235590174262844, 0.114508814744257199342353731044292225247093225};

// Machine independent implementation of pow(*,1./7.)
static double sqrt7(double a){
    double x = 1.;
    for (int k=0; k<20;k++){  // A smaller number should be ok too.
        double x6 = x*x*x*x*x*x;
        x += (a/x6-x)/7.;
    }
    return x;
}

static void free_dp7(struct reb_dp7* dp7){
    free(dp7->p0);
    free(dp7->p1);
    free(dp7->p2);
    free(dp7->p3);
    free(dp7->p4);
    free(dp7->p5);
    free(dp7->p6);
    dp7->p0 = NULL;
    dp7->p1 = NULL;
    dp7->p2 = NULL;
    dp7->p3 = NULL;
    dp7->p4 = NULL;
    dp7->p5 = NULL;
    dp7->p6 = NULL;
}
static void clear_dp7(struct reb_dp7* const dp7, const int N3){
    for (int k=0;k<N3;k++){
        dp7->p0[k] = 0.;
        dp7->p1[k] = 0.;
        dp7->p2[k] = 0.;
        dp7->p3[k] = 0.;
        dp7->p4[k] = 0.;
        dp7->p5[k] = 0.;
        dp7->p6[k] = 0.;
    }
}
static void realloc_dp7(struct reb_dp7* const dp7, const int N3){
    dp7->p0 = realloc(dp7->p0,sizeof(double)*N3);
    dp7->p1 = realloc(dp7->p1,sizeof(double)*N3);
    dp7->p2 = realloc(dp7->p2,sizeof(double)*N3);
    dp7->p3 = realloc(dp7->p3,sizeof(double)*N3);
    dp7->p4 = realloc(dp7->p4,sizeof(double)*N3);
    dp7->p5 = realloc(dp7->p5,sizeof(double)*N3);
    dp7->p6 = realloc(dp7->p6,sizeof(double)*N3);
    clear_dp7(dp7,N3);
}

static struct reb_dpconst7 dpcast(struct reb_dp7 dp){
    struct reb_dpconst7 dpc = {
        .p0 = dp.p0, 
        .p1 = dp.p1, 
        .p2 = dp.p2, 
        .p3 = dp.p3, 
        .p4 = dp.p4, 
        .p5 = dp.p5, 
        .p6 = dp.p6, 
    };
    return dpc;
}

static inline void add_cs(double* p, double* csp, double inp){
    const double y = inp - *csp;
    const double t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

void reb_integrator_ias15_alloc(struct reb_simulation* r){
    int N3;
    if (r->integrator==REB_INTEGRATOR_MERCURIUS){
        N3 = 3*r->ri_mercurius.encounterN;// mercurius close encounter
    }else{ 
        N3 = 3*r->N;
    }
    if (N3 > r->ri_ias15.allocatedN) {
        realloc_dp7(&(r->ri_ias15.g),N3);
        realloc_dp7(&(r->ri_ias15.b),N3);
        realloc_dp7(&(r->ri_ias15.csb),N3);
        realloc_dp7(&(r->ri_ias15.e),N3);
        realloc_dp7(&(r->ri_ias15.br),N3);
        realloc_dp7(&(r->ri_ias15.er),N3);
        r->ri_ias15.at = realloc(r->ri_ias15.at,sizeof(double)*N3);
        r->ri_ias15.x0 = realloc(r->ri_ias15.x0,sizeof(double)*N3);
        r->ri_ias15.v0 = realloc(r->ri_ias15.v0,sizeof(double)*N3);
        r->ri_ias15.a0 = realloc(r->ri_ias15.a0,sizeof(double)*N3);
        r->ri_ias15.csx= realloc(r->ri_ias15.csx,sizeof(double)*N3);
        r->ri_ias15.csv= realloc(r->ri_ias15.csv,sizeof(double)*N3);
        r->ri_ias15.csa0 = realloc(r->ri_ias15.csa0,sizeof(double)*N3);
        double* restrict const csx = r->ri_ias15.csx; 
        double* restrict const csv = r->ri_ias15.csv; 
        for (int i=0;i<N3;i++){
            // Kill compensated summation coefficients
            csx[i] = 0;
            csv[i] = 0;
        }
        r->ri_ias15.allocatedN = N3;
    }
    if (N3/3 > r->ri_ias15.map_allocated_N){
        r->ri_ias15.map = realloc(r->ri_ias15.map,sizeof(int)*(N3/3));
        for (int i=0;i<N3/3;i++){
            r->ri_ias15.map[i] = i;
        }
        r->ri_ias15.map_allocated_N = N3/3;
    }

}
 
// Does the actual timestep.
static int reb_integrator_ias15_step(struct reb_simulation* r) {
    reb_integrator_ias15_alloc(r);

    struct reb_particle* const particles = r->particles;
    int N;
    int* map; // this map allow for integrating only a selection of particles 
    if (r->integrator==REB_INTEGRATOR_MERCURIUS){// mercurius close encounter
        N = r->ri_mercurius.encounterN;
        map = r->ri_mercurius.encounter_map;
        if (map==NULL){
            reb_error(r, "Cannot access MERCURIUS map from IAS15.");
            return 0;
        }
    }else{ 
        N = r->N;
        map = r->ri_ias15.map; // identity map
    }
    const int N3 = 3*N;
    
    // reb_update_acceleration(); // Not needed. Forces are already calculated in main routine.
    
    // New order is better. Old order for backwards compatibility. 
    const unsigned int neworder = r->ri_ias15.neworder;  
    double s[9];                // Summation coefficients 
    double* restrict const csx = r->ri_ias15.csx; 
    double* restrict const csv = r->ri_ias15.csv; 
    double* restrict const csa0 = r->ri_ias15.csa0; 
    double* restrict const at = r->ri_ias15.at; 
    double* restrict const x0 = r->ri_ias15.x0; 
    double* restrict const v0 = r->ri_ias15.v0; 
    double* restrict const a0 = r->ri_ias15.a0; 
    struct reb_vec3d* gravity_cs = r->gravity_cs; 
    const struct reb_dpconst7 g  = dpcast(r->ri_ias15.g);
    const struct reb_dpconst7 e  = dpcast(r->ri_ias15.e);
    const struct reb_dpconst7 b  = dpcast(r->ri_ias15.b);
    const struct reb_dpconst7 csb= dpcast(r->ri_ias15.csb);
    const struct reb_dpconst7 er = dpcast(r->ri_ias15.er);
    const struct reb_dpconst7 br = dpcast(r->ri_ias15.br);
    for(int k=0;k<N;k++) {
        int mk = map[k];
        x0[3*k]   = particles[mk].x;
        x0[3*k+1] = particles[mk].y;
        x0[3*k+2] = particles[mk].z;
        v0[3*k]   = particles[mk].vx;
        v0[3*k+1] = particles[mk].vy;
        v0[3*k+2] = particles[mk].vz;
        a0[3*k]   = particles[mk].ax;
        a0[3*k+1] = particles[mk].ay; 
        a0[3*k+2] = particles[mk].az;
    }
    if (r->gravity==REB_GRAVITY_COMPENSATED){
        for(int k=0;k<N;k++) {
            int mk = map[k];
            csa0[3*k]   = gravity_cs[mk].x;
            csa0[3*k+1] = gravity_cs[mk].y;  
            csa0[3*k+2] = gravity_cs[mk].z;
        }
    }else{
        gravity_cs = (struct reb_vec3d*)csa0; // Always 0.
        for(int k=0;k<N3;k++) {
            csa0[k]   = 0;
        }
    }
    for (int k=0;k<N3;k++){
        // Memset might be faster!
        csb.p0[k] = 0.;
        csb.p1[k] = 0.;
        csb.p2[k] = 0.;
        csb.p3[k] = 0.;
        csb.p4[k] = 0.;
        csb.p5[k] = 0.;
        csb.p6[k] = 0.;
    }

    for(int k=0;k<N3;k++) {
        g.p0[k] = b.p6[k]*d[15] + b.p5[k]*d[10] + b.p4[k]*d[6] + b.p3[k]*d[3]  + b.p2[k]*d[1]  + b.p1[k]*d[0]  + b.p0[k];
        g.p1[k] = b.p6[k]*d[16] + b.p5[k]*d[11] + b.p4[k]*d[7] + b.p3[k]*d[4]  + b.p2[k]*d[2]  + b.p1[k];
        g.p2[k] = b.p6[k]*d[17] + b.p5[k]*d[12] + b.p4[k]*d[8] + b.p3[k]*d[5]  + b.p2[k];
        g.p3[k] = b.p6[k]*d[18] + b.p5[k]*d[13] + b.p4[k]*d[9] + b.p3[k];
        g.p4[k] = b.p6[k]*d[19] + b.p5[k]*d[14] + b.p4[k];
        g.p5[k] = b.p6[k]*d[20] + b.p5[k];
        g.p6[k] = b.p6[k];
    }

    double integrator_megno_thisdt = 0.;
    double integrator_megno_thisdt_init = 0.;
    if (r->calculate_megno){
        integrator_megno_thisdt_init = w[0]* r->t * reb_tools_megno_deltad_delta(r);
    }

    double t_beginning = r->t;
    double predictor_corrector_error = 1e300;
    double predictor_corrector_error_last = 2;
    int iterations = 0; 
    // Predictor corrector loop
    // Stops if one of the following conditions is satisfied: 
    //   1) predictor_corrector_error better than 1e-16 
    //   2) predictor_corrector_error starts to oscillate
    //   3) more than 12 iterations
    while(1){
        if(predictor_corrector_error<1e-16){
            break;
        }
        if(iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error){
            break;
        }
        if (iterations>=12){
            r->ri_ias15.iterations_max_exceeded++;
            const int integrator_iterations_warning = 10;
            if (r->ri_ias15.iterations_max_exceeded==integrator_iterations_warning ){
                reb_warning(r, "At least 10 predictor corrector loops in IAS15 did not converge. This is typically an indication of the timestep being too large.");
            }
            break;                              // Quit predictor corrector loop
        }
        predictor_corrector_error_last = predictor_corrector_error;
        predictor_corrector_error = 0;
        iterations++;

        integrator_megno_thisdt = integrator_megno_thisdt_init;

        for(int n=1;n<8;n++) {                          // Loop over interval using Gauss-Radau spacings

            if (!neworder){
                s[0] = r->dt * h[n];
                s[1] = s[0] * s[0] / 2.;
                s[2] = s[1] * h[n] / 3.;
                s[3] = s[2] * h[n] / 2.;
                s[4] = 3. * s[3] * h[n] / 5.;
                s[5] = 2. * s[4] * h[n] / 3.;
                s[6] = 5. * s[5] * h[n] / 7.;
                s[7] = 3. * s[6] * h[n] / 4.;
                s[8] = 7. * s[7] * h[n] / 9.;
            }
            
            r->t = t_beginning + r->dt * h[n];

            // Prepare particles arrays for force calculation
            for(int i=0;i<N;i++) {                      // Predict positions at interval n using b values
                int mi = map[i];
                const int k0 = 3*i+0;
                const int k1 = 3*i+1;
                const int k2 = 3*i+2;

                double xk0;
                double xk1;
                double xk2;
                if (neworder){
                    xk0 = -csx[k0] + ((((((((b.p6[k0]*7.*h[n]/9. + b.p5[k0])*3.*h[n]/4. + b.p4[k0])*5.*h[n]/7. + b.p3[k0])*2.*h[n]/3. + b.p2[k0])*3.*h[n]/5. + b.p1[k0])*h[n]/2. + b.p0[k0])*h[n]/3. + a0[k0])*r->dt*h[n]/2. + v0[k0])*r->dt*h[n];
                    xk1 = -csx[k1] + ((((((((b.p6[k1]*7.*h[n]/9. + b.p5[k1])*3.*h[n]/4. + b.p4[k1])*5.*h[n]/7. + b.p3[k1])*2.*h[n]/3. + b.p2[k1])*3.*h[n]/5. + b.p1[k1])*h[n]/2. + b.p0[k1])*h[n]/3. + a0[k1])*r->dt*h[n]/2. + v0[k1])*r->dt*h[n];
                    xk2 = -csx[k2] + ((((((((b.p6[k2]*7.*h[n]/9. + b.p5[k2])*3.*h[n]/4. + b.p4[k2])*5.*h[n]/7. + b.p3[k2])*2.*h[n]/3. + b.p2[k2])*3.*h[n]/5. + b.p1[k2])*h[n]/2. + b.p0[k2])*h[n]/3. + a0[k2])*r->dt*h[n]/2. + v0[k2])*r->dt*h[n];
                }else{
                    xk0 = -csx[k0] + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*a0[k0] + s[0]*v0[k0] );
                    xk1 = -csx[k1] + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*a0[k1] + s[0]*v0[k1] );
                    xk2 = -csx[k2] + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*a0[k2] + s[0]*v0[k2] );
                }
                particles[mi].x = xk0 + x0[k0];
                particles[mi].y = xk1 + x0[k1];
                particles[mi].z = xk2 + x0[k2];
            }
            if (r->calculate_megno || (r->additional_forces && r->force_is_velocity_dependent)){
                if (!neworder){
                    s[0] = r->dt * h[n];
                    s[1] =      s[0] * h[n] / 2.;
                    s[2] = 2. * s[1] * h[n] / 3.;
                    s[3] = 3. * s[2] * h[n] / 4.;
                    s[4] = 4. * s[3] * h[n] / 5.;
                    s[5] = 5. * s[4] * h[n] / 6.;
                    s[6] = 6. * s[5] * h[n] / 7.;
                    s[7] = 7. * s[6] * h[n] / 8.;
                }

                for(int i=0;i<N;i++) {                  // Predict velocities at interval n using b values
                    int mi = map[i];
                    const int k0 = 3*i+0;
                    const int k1 = 3*i+1;
                    const int k2 = 3*i+2;

                    double vk0;
                    double vk1;
                    double vk2;
                    if (neworder){
                        vk0 =  -csv[k0] + (((((((b.p6[k0]*7.*h[n]/8. + b.p5[k0])*6.*h[n]/7. + b.p4[k0])*5.*h[n]/6. + b.p3[k0])*4.*h[n]/5. + b.p2[k0])*3.*h[n]/4. + b.p1[k0])*2.*h[n]/3. + b.p0[k0])*h[n]/2. + a0[k0])*r->dt*h[n];
                        vk1 =  -csv[k1] + (((((((b.p6[k1]*7.*h[n]/8. + b.p5[k1])*6.*h[n]/7. + b.p4[k1])*5.*h[n]/6. + b.p3[k1])*4.*h[n]/5. + b.p2[k1])*3.*h[n]/4. + b.p1[k1])*2.*h[n]/3. + b.p0[k1])*h[n]/2. + a0[k1])*r->dt*h[n];
                        vk2 =  -csv[k2] + (((((((b.p6[k2]*7.*h[n]/8. + b.p5[k2])*6.*h[n]/7. + b.p4[k2])*5.*h[n]/6. + b.p3[k2])*4.*h[n]/5. + b.p2[k2])*3.*h[n]/4. + b.p1[k2])*2.*h[n]/3. + b.p0[k2])*h[n]/2. + a0[k2])*r->dt*h[n];
                    }else{
                        vk0 =  -csv[k0] + s[7]*b.p6[k0] + s[6]*b.p5[k0] + s[5]*b.p4[k0] + s[4]*b.p3[k0] + s[3]*b.p2[k0] + s[2]*b.p1[k0] + s[1]*b.p0[k0] + s[0]*a0[k0];
                        vk1 =  -csv[k1] + s[7]*b.p6[k1] + s[6]*b.p5[k1] + s[5]*b.p4[k1] + s[4]*b.p3[k1] + s[3]*b.p2[k1] + s[2]*b.p1[k1] + s[1]*b.p0[k1] + s[0]*a0[k1];
                        vk2 =  -csv[k2] + s[7]*b.p6[k2] + s[6]*b.p5[k2] + s[5]*b.p4[k2] + s[4]*b.p3[k2] + s[3]*b.p2[k2] + s[2]*b.p1[k2] + s[1]*b.p0[k2] + s[0]*a0[k2];
                    }
                    particles[mi].vx = vk0 + v0[k0];
                    particles[mi].vy = vk1 + v0[k1];
                    particles[mi].vz = vk2 + v0[k2];
                }
            }


            reb_update_acceleration(r);             // Calculate forces at interval n
            if (r->calculate_megno){
                integrator_megno_thisdt += w[n] * r->t * reb_tools_megno_deltad_delta(r);
            }

            for(int k=0;k<N;++k) {
                int mk = map[k];
                at[3*k]   = particles[mk].ax;
                at[3*k+1] = particles[mk].ay;  
                at[3*k+2] = particles[mk].az;
            }
            switch (n) {                            // Improve b and g values
                case 1: 
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p0[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p0[k]  = gk/rr[0];
                        add_cs(&(b.p0[k]), &(csb.p0[k]), g.p0[k]-tmp);
                    } break;
                case 2: 
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p1[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p1[k] = (gk/rr[1] - g.p0[k])/rr[2];
                        tmp = g.p1[k] - tmp;
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[0]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp);
                    } break;
                case 3: 
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p2[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p2[k] = ((gk/rr[3] - g.p0[k])/rr[4] - g.p1[k])/rr[5];
                        tmp = g.p2[k] - tmp;
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[1]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[2]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp);
                    } break;
                case 4:
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p3[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p3[k] = (((gk/rr[6] - g.p0[k])/rr[7] - g.p1[k])/rr[8] - g.p2[k])/rr[9];
                        tmp = g.p3[k] - tmp;
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[3]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[4]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[5]);
                        add_cs(&(b.p3[k]), &(csb.p3[k]), tmp);
                    } break;
                case 5:
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p4[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p4[k] = ((((gk/rr[10] - g.p0[k])/rr[11] - g.p1[k])/rr[12] - g.p2[k])/rr[13] - g.p3[k])/rr[14];
                        tmp = g.p4[k] - tmp;
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[6]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[7]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[8]);
                        add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[9]);
                        add_cs(&(b.p4[k]), &(csb.p4[k]), tmp);
                    } break;
                case 6:
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p5[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p5[k] = (((((gk/rr[15] - g.p0[k])/rr[16] - g.p1[k])/rr[17] - g.p2[k])/rr[18] - g.p3[k])/rr[19] - g.p4[k])/rr[20];
                        tmp = g.p5[k] - tmp;
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[10]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[11]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[12]);
                        add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[13]);
                        add_cs(&(b.p4[k]), &(csb.p4[k]), tmp * c[14]);
                        add_cs(&(b.p5[k]), &(csb.p5[k]), tmp);
                    } break;
                case 7:
                {
                    double maxak = 0.0;
                    double maxb6ktmp = 0.0;
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p6[k];
                        double gk = at[k];
                        double gk_cs = ((double*)(gravity_cs))[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p6[k] = ((((((gk/rr[21] - g.p0[k])/rr[22] - g.p1[k])/rr[23] - g.p2[k])/rr[24] - g.p3[k])/rr[25] - g.p4[k])/rr[26] - g.p5[k])/rr[27];
                        tmp = g.p6[k] - tmp;    
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[15]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[16]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[17]);
                        add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[18]);
                        add_cs(&(b.p4[k]), &(csb.p4[k]), tmp * c[19]);
                        add_cs(&(b.p5[k]), &(csb.p5[k]), tmp * c[20]);
                        add_cs(&(b.p6[k]), &(csb.p6[k]), tmp);
                        
                        // Monitor change in b.p6[k] relative to at[k]. The predictor corrector scheme is converged if it is close to 0.
                        if (r->ri_ias15.epsilon_global){
                            const double ak  = fabs(at[k]);
                            if (isnormal(ak) && ak>maxak){
                                maxak = ak;
                            }
                            const double b6ktmp = fabs(tmp);  // change of b6ktmp coefficient
                            if (isnormal(b6ktmp) && b6ktmp>maxb6ktmp){
                                maxb6ktmp = b6ktmp;
                            }
                        }else{
                            const double ak  = at[k];
                            const double b6ktmp = tmp; 
                            const double errork = fabs(b6ktmp/ak);
                            if (isnormal(errork) && errork>predictor_corrector_error){
                                predictor_corrector_error = errork;
                            }
                        }
                    } 
                    if (r->ri_ias15.epsilon_global){
                        predictor_corrector_error = maxb6ktmp/maxak;
                    }
                    
                    break;
                }
            }
        }
    }
    // Set time back to initial value (will be updated below) 
    r->t = t_beginning;
    // Find new timestep
    const double dt_done = r->dt;
    
    if (r->ri_ias15.epsilon>0){
        // Estimate error (given by last term in series expansion) 
        // There are two options:
        // r->ri_ias15.epsilon_global==1  (default)
        //   First, we determine the maximum acceleration and the maximum of the last term in the series. 
        //   Then, the two are divided.
        // r->ri_ias15.epsilon_global==0
        //   Here, the fractional error is calculated for each particle individually and we use the maximum of the fractional error.
        //   This might fail in cases where a particle does not experience any (physical) acceleration besides roundoff errors. 
        double integrator_error = 0.0;
        unsigned int Nreal = N - r->N_var;
        if (r->ri_ias15.neworder==0){
            Nreal = N; // backwards compatibility. will be removed in the future. 
        }
        if (r->ri_ias15.epsilon_global){
            double maxak = 0.0;
            double maxb6k = 0.0;
            for(int i=0;i<Nreal;i++){ // Looping over all particles and all 3 components of the acceleration. 
                // Note: Before December 2020, N-N_var, was simply N. This change should make timestep choices during
                // close encounters more stable if variational particles are present.
                int mi = map[i];
                const double v2 = particles[mi].vx*particles[mi].vx+particles[mi].vy*particles[mi].vy+particles[mi].vz*particles[mi].vz;
                const double x2 = particles[mi].x*particles[mi].x+particles[mi].y*particles[mi].y+particles[mi].z*particles[mi].z;
                // Skip slowly varying accelerations
                if (fabs(v2*r->dt*r->dt/x2) < 1e-16) continue;
                for(int k=3*i;k<3*(i+1);k++) { 
                    const double ak  = fabs(at[k]);
                    if (isnormal(ak) && ak>maxak){
                        maxak = ak;
                    }
                    const double b6k = fabs(b.p6[k]); 
                    if (isnormal(b6k) && b6k>maxb6k){
                        maxb6k = b6k;
                    }
                }
            }
            integrator_error = maxb6k/maxak;
        }else{
            for(int k=0;k<N3;k++) {
                const double ak  = at[k];
                const double b6k = b.p6[k]; 
                const double errork = fabs(b6k/ak);
                if (isnormal(errork) && errork>integrator_error){
                    integrator_error = errork;
                }
            }
        }

        double dt_new;
        if  (isnormal(integrator_error)){   
            // if error estimate is available increase by more educated guess
            dt_new = sqrt7(r->ri_ias15.epsilon/integrator_error)*dt_done;
        }else{                  // In the rare case that the error estimate doesn't give a finite number (e.g. when all forces accidentally cancel up to machine precission).
            dt_new = dt_done/safety_factor; // by default, increase timestep a little
        }
        
        if (fabs(dt_new)<r->ri_ias15.min_dt) dt_new = copysign(r->ri_ias15.min_dt,dt_new);
        
        if (fabs(dt_new/dt_done) < safety_factor) { // New timestep is significantly smaller.
            // Reset particles
            for(int k=0;k<N;++k) {
                int mk = map[k];
                particles[mk].x = x0[3*k+0]; // Set inital position
                particles[mk].y = x0[3*k+1];
                particles[mk].z = x0[3*k+2];

                particles[mk].vx = v0[3*k+0];    // Set inital velocity
                particles[mk].vy = v0[3*k+1];
                particles[mk].vz = v0[3*k+2];
                
                particles[mk].ax = a0[3*k+0];    // Set inital acceleration
                particles[mk].ay = a0[3*k+1];
                particles[mk].az = a0[3*k+2];
            }
            r->dt = dt_new;
            if (r->dt_last_done!=0.){       // Do not predict next e/b values if this is the first time step.
                double ratio = r->dt/r->dt_last_done;
                predict_next_step(ratio, N3, er, br, e, b);
            }
            
            return 0; // Step rejected. Do again. 
        }       
        if (fabs(dt_new/dt_done) > 1.0) {   // New timestep is larger.
            if (dt_new/dt_done > 1./safety_factor) dt_new = dt_done /safety_factor; // Don't increase the timestep by too much compared to the last one.
        }
        r->dt = dt_new;
    }

    // Find new position and velocity values at end of the sequence
    double dt_done2 = dt_done*dt_done;
    for(int k=0;k<N3;++k) {
        if (neworder){
            // Note: dt_done*dt_done is not precalculated to avoid 
            //       biased round-off errors when a fixed timestep is used.
            add_cs(&(x0[k]), &(csx[k]), b.p6[k]/72.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p5[k]/56.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p4[k]/42.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p3[k]/30.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p2[k]/20.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p1[k]/12.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p0[k]/6.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), a0[k]/2.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), v0[k]*dt_done);
        }else{
            add_cs(&(x0[k]), &(csx[k]), b.p6[k]/72.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p5[k]/56.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p4[k]/42.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p3[k]/30.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p2[k]/20.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p1[k]/12.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), b.p0[k]/6.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), a0[k]/2.*dt_done2);
            add_cs(&(x0[k]), &(csx[k]), v0[k]*dt_done);
        }
        {
            add_cs(&(v0[k]), &(csv[k]), b.p6[k]/8.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p5[k]/7.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p4[k]/6.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p3[k]/5.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p2[k]/4.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p1[k]/3.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p0[k]/2.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), a0[k]*dt_done);
        }
    }

    r->t += dt_done;
    r->dt_last_done = dt_done;

    if (r->calculate_megno){
        double dY = dt_done*integrator_megno_thisdt;
        reb_tools_megno_update(r, dY);
    }

    // Swap particle buffers
    for(int k=0;k<N;++k) {
        int mk = map[k];
        particles[mk].x = x0[3*k+0]; // Set final position
        particles[mk].y = x0[3*k+1];
        particles[mk].z = x0[3*k+2];

        particles[mk].vx = v0[3*k+0];    // Set final velocity
        particles[mk].vy = v0[3*k+1];
        particles[mk].vz = v0[3*k+2];
    }
    copybuffers(e,er,N3);       
    copybuffers(b,br,N3);       
    double ratio = r->dt/dt_done;
    predict_next_step(ratio, N3, e, b, e, b);
    return 1; // Success.
}

static void predict_next_step(double ratio, int N3,  const struct reb_dpconst7 _e, const struct reb_dpconst7 _b, const struct reb_dpconst7 e, const struct reb_dpconst7 b){
    if (ratio>20.){
        // Do not predict if stepsize increase is very large. 
        for(int k=0;k<N3;++k) {
            e.p0[k] = 0.; e.p1[k] = 0.; e.p2[k] = 0.; e.p3[k] = 0.; e.p4[k] = 0.; e.p5[k] = 0.; e.p6[k] = 0.;
            b.p0[k] = 0.; b.p1[k] = 0.; b.p2[k] = 0.; b.p3[k] = 0.; b.p4[k] = 0.; b.p5[k] = 0.; b.p6[k] = 0.;
        }
    }else{
        // Predict new B values to use at the start of the next sequence. The predicted
        // values from the last call are saved as E. The correction, BD, between the
        // actual and predicted values of B is applied in advance as a correction.
        //
        const double q1 = ratio;
        const double q2 = q1 * q1;
        const double q3 = q1 * q2;
        const double q4 = q2 * q2;
        const double q5 = q2 * q3;
        const double q6 = q3 * q3;
        const double q7 = q3 * q4;

        for(int k=0;k<N3;++k) {
            double be0 = _b.p0[k] - _e.p0[k];
            double be1 = _b.p1[k] - _e.p1[k];
            double be2 = _b.p2[k] - _e.p2[k];
            double be3 = _b.p3[k] - _e.p3[k];
            double be4 = _b.p4[k] - _e.p4[k];
            double be5 = _b.p5[k] - _e.p5[k];
            double be6 = _b.p6[k] - _e.p6[k];


            e.p0[k] = q1*(_b.p6[k]* 7.0 + _b.p5[k]* 6.0 + _b.p4[k]* 5.0 + _b.p3[k]* 4.0 + _b.p2[k]* 3.0 + _b.p1[k]*2.0 + _b.p0[k]);
            e.p1[k] = q2*(_b.p6[k]*21.0 + _b.p5[k]*15.0 + _b.p4[k]*10.0 + _b.p3[k]* 6.0 + _b.p2[k]* 3.0 + _b.p1[k]);
            e.p2[k] = q3*(_b.p6[k]*35.0 + _b.p5[k]*20.0 + _b.p4[k]*10.0 + _b.p3[k]* 4.0 + _b.p2[k]);
            e.p3[k] = q4*(_b.p6[k]*35.0 + _b.p5[k]*15.0 + _b.p4[k]* 5.0 + _b.p3[k]);
            e.p4[k] = q5*(_b.p6[k]*21.0 + _b.p5[k]* 6.0 + _b.p4[k]);
            e.p5[k] = q6*(_b.p6[k]* 7.0 + _b.p5[k]);
            e.p6[k] = q7* _b.p6[k];

            b.p0[k] = e.p0[k] + be0;
            b.p1[k] = e.p1[k] + be1;
            b.p2[k] = e.p2[k] + be2;
            b.p3[k] = e.p3[k] + be3;
            b.p4[k] = e.p4[k] + be4;
            b.p5[k] = e.p5[k] + be5;
            b.p6[k] = e.p6[k] + be6;
        }
    }
}

static void copybuffers(const struct reb_dpconst7 _a, const struct reb_dpconst7 _b, int N3){
    for (int i=0;i<N3;i++){ 
        _b.p0[i] = _a.p0[i];
        _b.p1[i] = _a.p1[i];
        _b.p2[i] = _a.p2[i];
        _b.p3[i] = _a.p3[i];
        _b.p4[i] = _a.p4[i];
        _b.p5[i] = _a.p5[i];
        _b.p6[i] = _a.p6[i];
    }
// The above code seems faster than the code below, probably due to some compiler optimizations. 
//  for (int i=0;i<7;i++){  
//      memcpy(_b[i],_a[i], sizeof(double)*N3);
//  }
}

// Do nothing here. This is only used in a leapfrog-like DKD integrator. IAS15 performs one complete timestep.
void reb_integrator_ias15_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
}

void reb_integrator_ias15_part2(struct reb_simulation* r){
#ifdef GENERATE_CONSTANTS
    integrator_generate_constants();
#endif  // GENERATE_CONSTANTS
    // Try until a step was successful.
    while(!reb_integrator_ias15_step(r));
}

void reb_integrator_ias15_synchronize(struct reb_simulation* r){
}
void reb_integrator_ias15_clear(struct reb_simulation* r){
    const int N3 = r->ri_ias15.allocatedN;
    if (N3){
        clear_dp7(&(r->ri_ias15.g),N3);
        clear_dp7(&(r->ri_ias15.e),N3);
        clear_dp7(&(r->ri_ias15.b),N3);
        clear_dp7(&(r->ri_ias15.csb),N3);
        clear_dp7(&(r->ri_ias15.er),N3);
        clear_dp7(&(r->ri_ias15.br),N3);
        
        double* restrict const csx = r->ri_ias15.csx; 
        double* restrict const csv = r->ri_ias15.csv; 
        for (int i=0;i<N3;i++){
            // Kill compensated summation coefficients
            csx[i] = 0;
            csv[i] = 0;
        }
    }
}

void reb_integrator_ias15_reset(struct reb_simulation* r){
    r->ri_ias15.allocatedN  = 0;
    r->ri_ias15.map_allocated_N  = 0;
    free_dp7(&(r->ri_ias15.g));
    free_dp7(&(r->ri_ias15.e));
    free_dp7(&(r->ri_ias15.b));
    free_dp7(&(r->ri_ias15.csb));
    free_dp7(&(r->ri_ias15.er));
    free_dp7(&(r->ri_ias15.br));
    free(r->ri_ias15.at);
    r->ri_ias15.at =  NULL;
    free(r->ri_ias15.x0);
    r->ri_ias15.x0 =  NULL;
    free(r->ri_ias15.v0);
    r->ri_ias15.v0 =  NULL;
    free(r->ri_ias15.a0);
    r->ri_ias15.a0 =  NULL;
    free(r->ri_ias15.csx);
    r->ri_ias15.csx=  NULL;
    free(r->ri_ias15.csv);
    r->ri_ias15.csv=  NULL;
    free(r->ri_ias15.csa0);
    r->ri_ias15.csa0 =  NULL;
    free(r->ri_ias15.map);
    r->ri_ias15.map =  NULL;
}

#ifdef GENERATE_CONSTANTS
void integrator_generate_constants(void){
    printf("Generaring constants.\n\n");
    mpf_set_default_prec(512);
    mpf_t* _h = malloc(sizeof(mpf_t)*8);
    for (int i=0;i<8;i++){
        mpf_init(_h[i]);
    }
    mpf_t* _r = malloc(sizeof(mpf_t)*28);
    for (int i=0;i<28;i++){
        mpf_init(_r[i]);
    }
    mpf_t* _c = malloc(sizeof(mpf_t)*21);
    mpf_t* _d = malloc(sizeof(mpf_t)*21);
    for (int i=0;i<21;i++){
        mpf_init(_c[i]);
        mpf_init(_d[i]);
    }
    mpf_set_str(_h[0],"0.0",10);
    mpf_set_str(_h[1],"0.0562625605369221464656521910318",10);
    mpf_set_str(_h[2],"0.180240691736892364987579942780",10);
    mpf_set_str(_h[3],"0.352624717113169637373907769648",10);
    mpf_set_str(_h[4],"0.547153626330555383001448554766",10);
    mpf_set_str(_h[5],"0.734210177215410531523210605558", 10);
    mpf_set_str(_h[6],"0.885320946839095768090359771030",10);
    mpf_set_str(_h[7],"0.977520613561287501891174488626",10);

    int l=0;
    for (int j=1;j<8;++j) {
        for(int k=0;k<j;++k) {
            // rr[l] = h[j] - h[k];
            mpf_sub(_r[l],_h[j],_h[k]);
            ++l;
        }
    }
    //c[0] = -h[1];
    mpf_neg(_c[0],_h[1]);
    //d[0] =  h[1];
    mpf_set(_d[0],_h[1]);
    l=0;
    for (int j=2;j<7;++j) {
        ++l;
        // c[l] = -h[j] * c[l-j+1];
        mpf_mul(_c[l], _h[j], _c[l-j+1]);
        mpf_neg(_c[l], _c[l]);
        //d[l] =  h[1] * d[l-j+1];
        mpf_mul(_d[l], _h[1], _d[l-j+1]);
        for(int k=2;k<j;++k) {
            ++l;
            //c[l] = c[l-j] - h[j] * c[l-j+1];
            mpf_mul(_c[l], _h[j], _c[l-j+1]);
            mpf_sub(_c[l], _c[l-j], _c[l]);
            //d[l] = d[l-j] + h[k] * d[l-j+1];
            mpf_mul(_d[l], _h[k], _d[l-j+1]);
            mpf_add(_d[l], _d[l-j], _d[l]);
        }
        ++l;
        //c[l] = c[l-j] - h[j];
        mpf_sub(_c[l], _c[l-j], _h[j]);
        //d[l] = d[l-j] + h[j]; 
        mpf_add(_d[l], _d[l-j], _h[j]);
    }

    // Output   
    printf("double rr[28] = {");
    for (int i=0;i<28;i++){
         gmp_printf ("%.*Ff", 25, _r[i]);
         if (i!=27) printf(", ");
    }
    printf("};\n");
    printf("double c[21] = {");
    for (int i=0;i<21;i++){
         gmp_printf ("%.*Ff", 25, _c[i]);
         if (i!=20) printf(", ");
    }
    printf("};\n");
    printf("double d[21] = {");
    for (int i=0;i<21;i++){
         gmp_printf ("%.*Ff", 25, _d[i]);
         if (i!=20) printf(", ");
    }
    printf("};\n");
    exit(0);
}
#endif // GENERATE_CONSTANTS
