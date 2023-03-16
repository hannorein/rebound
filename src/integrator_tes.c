/**
 * @file 	integrator_tes.c
 * @brief 	Terrestrial Exoplanet Simulator (TES) integration scheme
 * @author 	Peter Bartram <p.bartram@soton.ac.uk>
 * 
 * @section 	LICENSE
 * Copyright (c) 2022  <Peter Bartram, Alexander Wittig>
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
#include <stdlib.h>
#include <string.h>
#include "rebound.h"
#include "integrator_tes.h"

#define stages 7
#define FINAL_STAGE_INDEX 8
#define MAX_STEP_SIZE_GROWTH 4.0
#define MIN_STEP_SIZE 1E-5
#define OSCULATING_ORBIT_SLOTS 9  // stages + 2 for t=0 and t=1.
#define MAX_ITERATIONS 12
#define MIN_ITERATIONS 2
#define PI 3.141592653589793238462643383279
#define PI_SQ_X4 (4.0*PI*PI)
#define MAX_NEWTON_ITERATIONS 50
#define STUMPF_ITERATIONS 13 

// Coefficients
static const double hArr[9] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1.0};
static const double rr_inv[28] = {17.7738089140780008407526624, 5.5481367185372165056928203, 8.0659386483818866885371223, 2.8358760786444386782520107, 3.3742499769626352599420361, 5.8010015592640614823286804, 1.8276402675175978297946079, 2.0371118353585847827949161, 2.7254422118082262837742729, 5.1406241058109342286363203, 1.3620078160624694969370006, 1.4750402175604115479218482, 1.8051535801402512604391149, 2.6206449263870350811541814, 5.3459768998711075141214895, 1.1295338753367899027322862, 1.2061876660584456166252037, 1.4182782637347391537713785, 1.8772424961868100972169920, 2.9571160172904557478071039, 6.6176620137024244874471300, 1.0229963298234867458386119, 1.0854721939386423840467243, 1.2542646222818777659905423, 1.6002665494908162609916716, 2.3235983002196942228325344, 4.1099757783445590862385765, 10.8460261902368446847064289};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};
static const double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588};
static double invfactorial[35] = {1.0,1.0,0.5,0.166666666666666666666666666666667,0.0416666666666666666666666666666667,
                            0.00833333333333333333333333333333333,0.00138888888888888888888888888888889,0.000198412698412698412698412698412698,
                            0.0000248015873015873015873015873015873,0.00000275573192239858906525573192239859,0.000000275573192239858906525573192239859,
                            0.0000000250521083854417187750521083854417,0.00000000208767569878680989792100903212014,0.000000000160590438368216145993923771701549,
                            1.14707455977297247138516979786821e-11,7.64716373181981647590113198578807e-13,4.77947733238738529743820749111754e-14,
                            2.81145725434552076319894558301032e-15,1.56192069685862264622163643500573e-16,8.22063524662432971695598123687228e-18,
                            4.11031762331216485847799061843614e-19,1.95729410633912612308475743735054e-20,8.89679139245057328674889744250247e-22,
                            3.86817017063068403771691193152281e-23,1.61173757109611834904871330480117e-24,6.44695028438447339619485321920469e-26,
                            2.47959626322479746007494354584796e-27,9.18368986379554614842571683647391e-29,3.27988923706983791015204172731211e-30,
                            1.13099628864477169315587645769383e-31,3.76998762881590564385292152564611e-33,1.21612504155351794962997468569229e-34,
                            3.80039075485474359259367089278841e-36,1.15163356207719502805868814932982e-37,3.38715753552116184723143573332301e-39};

typedef struct _controlVars_const {
    const double* const __restrict__ p0;
    const double* const __restrict__ p1;
    const double* const __restrict__ p2;
    const double* const __restrict__ p3;
    const double* const __restrict__ p4;
    const double* const __restrict__ p5;
    const double* const __restrict__ p6;
    uint32_t size;
}controlVars_const;

static controlVars_const control_vars_cast(controlVars in)
{
    controlVars_const out = {
        .p0 = in.p0, 
        .p1 = in.p1, 
        .p2 = in.p2, 
        .p3 = in.p3, 
        .p4 = in.p4, 
        .p5 = in.p5, 
        .p6 = in.p6, 
    };    
    return out;
}

// Macros for performance.
#define NORM(x, i) sqrt(x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2])
#define DOT(x, y, i) x[3*i+0]*y[3*i+0]+x[3*i+1]*y[3*i+1]+x[3*i+2]*y[3*i+2]
#define CROSS(a,b,c,i) \
	(a)[3*i+0] = (b)[3*i+1] * (c)[3*i+2] - (c)[3*i+1] * (b)[3*i+2]; \
	(a)[3*i+1] = (b)[3*i+2] * (c)[3*i+0] - (c)[3*i+2] * (b)[3*i+0]; \
	(a)[3*i+2] = (b)[3*i+0] * (c)[3*i+1] - (c)[3*i+0] * (b)[3*i+1];
#define CROSS_LOCAL(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


// Top level functions
static void reb_tes_init(struct reb_simulation* r, uint32_t z_n);
static void reb_tes_free(struct reb_simulation* r);

// Gravity functions
static void reb_dhem_perform_summation(struct reb_simulation* r, double * Q, double * P,
                               double * dQ, double * dP, uint32_t stageNumber);
static void reb_dhem_calc_osc_orbits_for_all_stages(struct reb_simulation* r, double t0, double h, const double * hArr, uint32_t z_stagesPerStep, uint32_t z_rebasis);
static void reb_dhem_init_osc_orbits(struct reb_simulation* r, double * Q, double * P, double t);
static void reb_dhem_rhs(struct reb_simulation* r, double const * __restrict__ const dQ, double const * __restrict__ const dP, double * __restrict__ const dQ_dot,
              double * __restrict__ const dP_dot, double * __restrict__ const dQ_ddot, double * __restrict__ const dP_ddot);
static void reb_dhem_rhs_wrapped(struct reb_simulation* r, double * dQ, double * dP, double * dQ_dot,
                      double * dP_dot, double * dQ_ddot, double * dP_ddot, uint32_t stageNumber,
                      double * cs1, double * cs2);      
static void reb_dhem_calc_osc_orbit_derivs(struct reb_simulation* r, double const * const __restrict__ Qosc, double const * const __restrict__ Posc, 
                                                      double * const __restrict__ Qosc_dot, double * const __restrict__ Posc_dot);
static uint32_t reb_dhem_rectify(struct reb_simulation* r, double t, double * Q, double * P,
                            double * dQ, double * dP, uint32_t * rectifiedArray, uint32_t stageNumber);
static double reb_find_min_particle_period(struct reb_simulation* r);
static void reb_dhem_init(struct reb_simulation* r, double z_rectificationPeriodDefault, uint32_t z_stagesPerStep);
static void reb_dhem_free(struct reb_simulation* r);

// Radau functions
static void reb_radau_init(struct reb_simulation* r);
static void reb_radau_free(struct reb_simulation* r);
static double reb_calc_stepsize(struct reb_simulation* r, double h, double hLast, double t);
static void reb_clear_rectified_b_fields(struct reb_simulation* r, controlVars * B, uint32_t * rectifiedArray);
static double reb_single_step(struct reb_simulation* r, double z_t, double dt, double dt_last_done);
static void reb_init_radau_step(struct reb_simulation* r);
static void reb_free_radau_step(struct reb_simulation* r);
static void reb_free_controlvars(controlVars * var);
static void reb_init_controlvars(controlVars * var, uint32_t size);
static void reb_clear_controlvars(controlVars * var);
static void reb_calc_predictors(double h, double hSample, double const * __restrict__ z_state0, double const * __restrict__ z_dState, 
                         double const * __restrict__ z_ddState, controlVars const * z_B, double * __restrict__ z_predictors, 
                         double const * __restrict__ z_csState, uint32_t const z_start, uint32_t const z_end);
static void reb_calc_predictors_1st_order(double h, double hSample, double * z_state0, double * z_dState, 
                                  controlVars * z_B, double * z_predictors, double *, uint32_t start, uint32_t length);

static void reb_update_state(double h, double * z_dState, double * z_ddState, controlVars * z_B, double * z_state0, double * z_csState, uint32_t z_start, uint32_t z_end);
static void reb_update_state_1st_order(double h, double * z_dState, controlVars * z_B, double * z_state0, double * z_csState, uint32_t z_start, uint32_t z_end);
static void reb_analytical_continuation(struct reb_simulation* r, controlVars * z_B, controlVars * z_Blast, const double h, const double h_new, const uint32_t * const rectificationArray);
static void reb_calc_g_from_b(struct reb_simulation* r);\
static void reb_calc_g_from_b_internal(controlVars * z_G, controlVars * z_B, uint32_t z_start, uint32_t z_end);
static void reb_radau_step(struct reb_simulation* r, uint32_t * z_fCalls, double t, double h);
static double reb_calc_step_error(struct reb_simulation* r, double h, double t);

// Universal variables functions
static void reb_init_uvars(struct reb_simulation* const r);
static void reb_free_uvars(struct reb_simulation* r);
static void reb_rebasis_osc_orbits(struct reb_simulation* r, double * z_Q, double * z_P, double z_t, uint32_t i);
static void reb_calc_osc_orbits(struct reb_simulation* r, double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis);
static void reb_apply_osc_orbit_corrector(struct reb_simulation* r, double **Xosc_map, double t, uint32_t z_stagesPerStep);                                            
static void reb_c_stumpff(double * cs_in, double z_in);
static double reb_solve_for_universal_anomaly(struct reb_simulation* const r, double dt, double h, uint32_t i, double * C);
static double reb_calc_osc_update_value(struct reb_simulation* const r, double X, double dt, uint32_t i, double * C);


static inline void add_cs(double* out, double* cs, double inp);




void reb_integrator_tes_part1(struct reb_simulation* r){
    // unused
    if (r->ri_tes.warnings==0){
        r->ri_tes.warnings++;
        reb_warning(r,"The TES integrator is a recent addition to REBOUND. It may contain bugs or produce unphysical results. If you encounter a problem, please open an issue on GitHub.");
    }

} 

void reb_integrator_tes_part2(struct reb_simulation* r){
    uint32_t N = r->N;

    // If the number of bodies or the mass of the star has changed then do a "hard" reset.
    if(r->ri_tes.allocated_N != N || r->ri_tes.mStar_last != r->particles[0].m)
    {
        r->ri_tes.mStar_last = r->particles[0].m;
        r->ri_tes.allocated_N = N;
        struct reb_particle* const particles = r->particles;
        
        reb_tes_init(r, N);

        // Convert from inertial to dh coords.
        reb_transformations_inertial_to_democraticheliocentric_posvel(particles, r->ri_tes.particles_dh, r->N, r->N);

        // Store the COM and it's velocity as they are treated differently by TES.
        r->ri_tes.COM[0] = r->ri_tes.particles_dh[0].x;
        r->ri_tes.COM[1] = r->ri_tes.particles_dh[0].y;
        r->ri_tes.COM[2] = r->ri_tes.particles_dh[0].z;

        r->ri_tes.COM_dot[0] = r->ri_tes.particles_dh[0].vx;
        r->ri_tes.COM_dot[1] = r->ri_tes.particles_dh[0].vy;
        r->ri_tes.COM_dot[2] = r->ri_tes.particles_dh[0].vz;        
        
        for(uint32_t i=1;i<N;i++) 
        {
            r->ri_tes.mass[i] =     r->ri_tes.particles_dh[i].m;
            r->ri_tes.Q_dh[3*i] =   r->ri_tes.particles_dh[i].x;
            r->ri_tes.Q_dh[3*i+1] = r->ri_tes.particles_dh[i].y;
            r->ri_tes.Q_dh[3*i+2] = r->ri_tes.particles_dh[i].z;
            r->ri_tes.P_dh[3*i] =   r->ri_tes.particles_dh[i].vx*r->ri_tes.particles_dh[i].m;
            r->ri_tes.P_dh[3*i+1] = r->ri_tes.particles_dh[i].vy*r->ri_tes.particles_dh[i].m;
            r->ri_tes.P_dh[3*i+2] = r->ri_tes.particles_dh[i].vz*r->ri_tes.particles_dh[i].m; 
        }
        r->ri_tes.mass[0] = particles[0].m; // Keep mass[0] as stellar mass.

        reb_init_uvars(r);
        reb_dhem_init(r, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
        reb_dhem_init_osc_orbits(r, r->ri_tes.Q_dh, r->ri_tes.P_dh, r->t);
        reb_radau_init(r);  

        // No step-size specified, calculate timestep from shortest orbital period object in the system.
        if(r->dt == 0.001 && r->dt_last_done == 0)
        {
            // Default to 100th of minimum period for the initial step size.                                  
            r->dt = reb_find_min_particle_period(r)/100.0;
        }   

        // No orbital period specified, get this data from osculating elements.
        if(r->ri_tes.orbital_period == 0.0)
        {            
            r->ri_tes.orbital_period = reb_find_min_particle_period(r);

            // Update the rectification times again due to period change.
            for(int32_t i = 0; i < r->N; i++)
            {
              r->ri_tes.rhs->rectifyTimeArray[i] = r->t + r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit;
              r->ri_tes.rhs->rectificationPeriod[i] = r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit;
            }
        }     
    }
    
    double dt_new = reb_single_step(r, r->t, r->dt, r->dt_last_done);

    // update timestep
    r->t+=r->dt;
    r->dt_last_done = r->dt;
    r->dt = dt_new;
}


void reb_integrator_tes_synchronize(struct reb_simulation* r){
    uint32_t N = r->N;

    if(r->ri_tes.allocated_N == N)
    {    
        reb_dhem_perform_summation(r, r->ri_tes.radau->Qout, r->ri_tes.radau->Pout, 
                                        r->ri_tes.radau->dQ, r->ri_tes.radau->dP, 8);
                    
        double * Q_out = r->ri_tes.radau->Qout;
        double * P_out = r->ri_tes.radau->Pout;
        double * m = r->ri_tes.mass;

        for(uint32_t i=1; i < N; i++)
        {
            r->ri_tes.particles_dh[i].x = Q_out[3*i];
            r->ri_tes.particles_dh[i].y = Q_out[3*i+1];
            r->ri_tes.particles_dh[i].z = Q_out[3*i+2];
            r->ri_tes.particles_dh[i].vx = P_out[3*i]/m[i];    
            r->ri_tes.particles_dh[i].vy = P_out[3*i+1]/m[i];
            r->ri_tes.particles_dh[i].vz = P_out[3*i+2]/m[i];

            r->ri_tes.particles_dh[i].m = m[i];
        }       

        // Update the output with the COM including drift.
        r->ri_tes.particles_dh[0].x = r->ri_tes.COM[0];
        r->ri_tes.particles_dh[0].y = r->ri_tes.COM[1];
        r->ri_tes.particles_dh[0].z = r->ri_tes.COM[2];

        r->ri_tes.particles_dh[0].vx = r->ri_tes.COM_dot[0];
        r->ri_tes.particles_dh[0].vy = r->ri_tes.COM_dot[1];
        r->ri_tes.particles_dh[0].vz = r->ri_tes.COM_dot[2];  

        r->ri_tes.particles_dh[0].m = r->ri_tes.rhs->mTotal;

        reb_transformations_democraticheliocentric_to_inertial_posvel(r->particles, r->ri_tes.particles_dh, r->N, r->N); 
    }         
}

void reb_integrator_tes_reset(struct reb_simulation* r){
    if(r->ri_tes.allocated_N != 0)
    {
      reb_free_uvars(r);
      reb_dhem_free(r);
      reb_radau_free(r);
      reb_tes_free(r);    
    }

}

void reb_integrator_tes_allocate_memory(struct reb_simulation* r)
{
    reb_tes_init(r, r->N);
    reb_init_uvars(r);
    reb_dhem_init(r, r->ri_tes.orbital_period/r->ri_tes.recti_per_orbit, 9);
    reb_radau_init(r);          
}


///////////////////////////////////////////////////////////////////////////////////
// Original TES file: Simulation.c
///////////////////////////////////////////////////////////////////////////////////
static void reb_tes_init(struct reb_simulation* r, uint32_t z_n)
{
  if(r->ri_tes.allocated_N != 0)
  {
    reb_tes_free(r);
  }
  
  // Set control variables to initial values
  r->ri_tes.stateVectorLength = 2*3*r->N;
  r->ri_tes.stateVectorSize = r->ri_tes.stateVectorLength * sizeof(double);
  r->ri_tes.controlVectorSize = r->N * sizeof(double);

  // Allocate memory
  r->ri_tes.mass = (double *)malloc(r->ri_tes.controlVectorSize);
  r->ri_tes.X_dh = (double *)malloc(r->ri_tes.stateVectorSize);
  r->ri_tes.particles_dh = (struct reb_particle*)malloc(sizeof(struct reb_particle)*r->N);
  r->ri_tes.Q_dh = r->ri_tes.X_dh;
  r->ri_tes.P_dh = &r->ri_tes.X_dh[r->ri_tes.stateVectorLength/2];

  // Ensure we are clean for each integration.
  memset(r->ri_tes.mass, 0, r->ri_tes.controlVectorSize);
  memset(r->ri_tes.X_dh, 0, r->ri_tes.stateVectorSize);
}


static void reb_tes_free(struct reb_simulation* r)
{
  free(r->ri_tes.mass);
  free(r->ri_tes.X_dh);
  free(r->ri_tes.particles_dh);
}


///////////////////////////////////////////////////////////////////////////////////
// Original TES file: dhem.c
///////////////////////////////////////////////////////////////////////////////////
static void reb_dhem_rhs_wrapped(struct reb_simulation* r, double * dQ, double * dP, double * dQ_dot,
              double * dP_dot, double * dQ_ddot, double * dP_ddot, uint32_t stageNumber,
              double * cs1, double * cs2)
{
  DHEM * dhem = r->ri_tes.rhs;
  // Set up the pointer to the previously calculated osculating orbit values.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*r->N];
  dhem->Xosc_dot = dhem->Xosc_dotArr[stageNumber];
  dhem->Qosc_dot = dhem->Xosc_dot;
  dhem->Posc_dot = &dhem->Qosc_dot[3*r->N];

  // cs vars
  dhem->Xosc_cs = dhem->XoscArr_cs[stageNumber];
  dhem->Qosc_cs = dhem->Xosc_cs;
  dhem->Posc_cs = &dhem->Qosc_cs[3*r->N];

  reb_dhem_rhs(r, dQ, dP, dQ_dot, dP_dot, dQ_ddot, dP_ddot);

}

static void reb_dhem_rhs(struct reb_simulation* r, double const * __restrict__ const dQ, double const * __restrict__ const dP, double * __restrict__ const dQ_dot,
              double * __restrict__ const dP_dot, double * __restrict__ const dQ_ddot, double * __restrict__ const dP_ddot)
{
  DHEM * dhem = r->ri_tes.rhs;
  // Not necessary but makes code more reable.
  const double * const __restrict__ m = dhem->m;
  const uint32_t n = r->N;
  const uint32_t n3 = 3*r->N;
  const double G = r->G;
  const double * const __restrict__ Qosc = dhem->Qosc;
  const double * const __restrict__ Posc = dhem->Posc;
  const double * const __restrict__ Posc_dot = dhem->Posc_dot;
  double * __restrict__ Q = dhem->Q;
  double * __restrict__ P = dhem->P;

  memset(dQ_ddot, 0, r->ri_tes.stateVectorSize/2);

  #pragma GCC ivdep
  for(uint32_t i = 3; i < n3; i++)
  {
      Q[i] = Qosc[i] + dQ[i];
      P[i] = Posc[i] + dP[i]; 
  }    

  double vCentral[3] = {0,0,0};

  #pragma GCC ivdep
  for(uint32_t i = 1; i < n; i++)
  {
    vCentral[0] += P[3*i];
    vCentral[1] += P[3*i+1];
    vCentral[2] += P[3*i+2];
  }

  #pragma GCC ivdep
  for(uint32_t i = 0; i < 3; i ++)
  {
    vCentral[i] *= dhem->m_inv[0];
  }

  #pragma GCC ivdep
  for(uint32_t i = 1; i < n; i++)
  {
      dQ_dot[3*i] =   (dP[3*i]   / m[i]) + vCentral[0];
      dQ_dot[3*i+1] = (dP[3*i+1] / m[i]) + vCentral[1];
      dQ_dot[3*i+2] = (dP[3*i+2] / m[i]) + vCentral[2];
  }

  #pragma GCC ivdep
  for(uint32_t i = 1; i < n; i++)
  {
    const double GMM = G*m[0]*m[i];
    
    const double dQx = dQ[3*i+0];
    const double dQy = dQ[3*i+1];
    const double dQz = dQ[3*i+2];
    
    const double Qx = Q[3*i+0];
    const double Qy = Q[3*i+1];
    const double Qz = Q[3*i+2];      
    
    // Calculate the osculating orbit contribution term.
    const double Qosc_norm = NORM(Qosc, i);
    const double Qosc_norm3 = Qosc_norm*Qosc_norm*Qosc_norm;
    const double GMM_Qosc_norm3Inv = GMM/Qosc_norm3;

    const double drx = dQx-2*Qx;
    const double dry = dQy-2*Qy;
    const double drz = dQz-2*Qz;
    const double q = (dQx*drx+dQy*dry+dQz*drz) /  (Qx*Qx+Qy*Qy+Qz*Qz);
    const double q1 = (1+q);
    const double q3 = q1*q1*q1;
    const double fq = -q*(3+3*q+q*q) / (1+sqrt(q3));
    const double GMM_Qosc_norm3Inv_fq = GMM_Qosc_norm3Inv*fq;

    dP_dot[3*i]   = -dQx*GMM_Qosc_norm3Inv + GMM_Qosc_norm3Inv_fq*Qx;
    dP_dot[3*i+1] = -dQy*GMM_Qosc_norm3Inv + GMM_Qosc_norm3Inv_fq*Qy;
    dP_dot[3*i+2] = -dQz*GMM_Qosc_norm3Inv + GMM_Qosc_norm3Inv_fq*Qz;
  }
  
  const double istart = 1;
  for(uint32_t i = istart; i < n; i++)
  {
      const double GM = G*m[i];
      #pragma GCC ivdep
      for(uint32_t j = istart; j < i; j++)
      {
          const double GMM = GM*m[j];
          const double dx = Q[3*j+0] - Q[3*i+0];
          const double dy = Q[3*j+1] - Q[3*i+1];
          const double dz = Q[3*j+2] - Q[3*i+2];

          const double sepNorm = sqrt(dx*dx+dy*dy+dz*dz);
          const double sepNorm2 = sepNorm*sepNorm;
          const double sepNorm3 = sepNorm2*sepNorm;

          const double GMM_SepNorm3Inv = (GMM/sepNorm3);

          dP_dot[3*i+0] += dx*GMM_SepNorm3Inv;
          dP_dot[3*i+1] += dy*GMM_SepNorm3Inv;
          dP_dot[3*i+2] += dz*GMM_SepNorm3Inv;

          dP_dot[3*j+0] -= dx*GMM_SepNorm3Inv;
          dP_dot[3*j+1] -= dy*GMM_SepNorm3Inv;
          dP_dot[3*j+2] -= dz*GMM_SepNorm3Inv;
      }
  }


  double vCentralDot[3] = {0,0,0};

  /** Calculate the momentum of the central body. **/
  #pragma GCC ivdep
  for(uint32_t i = 1; i < n; i++)
  {
    vCentralDot[0] += (Posc_dot[3*i] + dP_dot[3*i]);
    vCentralDot[1] += (Posc_dot[3*i+1] + dP_dot[3*i+1]);
    vCentralDot[2] += (Posc_dot[3*i+2] + dP_dot[3*i+2]);
  }

  #pragma GCC ivdep
  for(uint32_t i = 0; i < 3; i ++)
  {
    vCentralDot[i] *= dhem->m_inv[0];
  }

  #pragma GCC ivdep
  for(uint32_t i = 1; i < n; i++)
  {
    dQ_ddot[3*i]   += vCentralDot[0];
    dQ_ddot[3*i+1] += vCentralDot[1];
    dQ_ddot[3*i+2] += vCentralDot[2];
    
    dQ_ddot[3*i] += dP_dot[3*i] * dhem->m_inv[i];
    dQ_ddot[3*i+1] += dP_dot[3*i+1] * dhem->m_inv[i];
    dQ_ddot[3*i+2] += dP_dot[3*i+2] * dhem->m_inv[i];
  }
}

static void reb_dhem_perform_summation(struct reb_simulation* r, double * Q, double * P,
                               double * dQ, double * dP, uint32_t stageNumber)
{
  DHEM * dhem = r->ri_tes.rhs;
  // Need to setup a pointer to the osculating orbits in memory.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*r->N];

  // Always zero in our frame of reference.
  for(uint32_t i = 0; i < 3; i++)
  {
    Q[i] = 0;
    P[i] = 0;
  }

  for(int32_t i = 1; i < r->N; i++)
  {
    Q[3*i+0] = dhem->Qosc[3*i+0] + (dQ[3*i+0] + (dhem->Qosc_cs[3*i+0] + r->ri_tes.radau->cs_dq[3*i+0]));
    Q[3*i+1] = dhem->Qosc[3*i+1] + (dQ[3*i+1] + (dhem->Qosc_cs[3*i+1] + r->ri_tes.radau->cs_dq[3*i+1]));
    Q[3*i+2] = dhem->Qosc[3*i+2] + (dQ[3*i+2] + (dhem->Qosc_cs[3*i+2] + r->ri_tes.radau->cs_dq[3*i+2]));

    P[3*i+0] = dhem->Posc[3*i+0] + (dP[3*i+0] + (dhem->Posc_cs[3*i+0] + r->ri_tes.radau->cs_dp[3*i+0]));
    P[3*i+1] = dhem->Posc[3*i+1] + (dP[3*i+1] + (dhem->Posc_cs[3*i+1] + r->ri_tes.radau->cs_dp[3*i+1]));
    P[3*i+2] = dhem->Posc[3*i+2] + (dP[3*i+2] + (dhem->Posc_cs[3*i+2] + r->ri_tes.radau->cs_dp[3*i+2]));
  }

}

static void reb_dhem_init_osc_orbits(struct reb_simulation* r, double * Q, double * P, double t)
{
  for(int32_t i = 1; i < r->N; i++)
  {
    reb_rebasis_osc_orbits(r, Q, P, t, i);
  }
}

uint32_t reb_dhem_rectify(struct reb_simulation* r, double t, double * Q, double * P,
                            double * dQ, double * dP, uint32_t * rectifiedArray, uint32_t stageNumber)
{
  DHEM * dhem = r->ri_tes.rhs;
  uint32_t rectifyFlag = 0;
  uint32_t rectifiedCount = 0;
  double dQ_norm = 0;

  // Need to setup a pointer to the osculating orbits in memory.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*r->N];
  dhem->Xosc_dot = dhem->Xosc_dotArr[stageNumber];
  dhem->Qosc_dot = dhem->Xosc_dot;
  dhem->Posc_dot = &dhem->Qosc_dot[3*r->N];

  // CS variables
  dhem->Xosc_cs = dhem->XoscArr_cs[stageNumber];
  dhem->Qosc_cs = dhem->Xosc_cs;
  dhem->Posc_cs = &dhem->Qosc_cs[3*r->N];

  for(int32_t i = 1; i < r->N; i++)
  {
    rectifiedArray[3*i] = 0;
    rectifiedArray[3*i+1] = 0;
    rectifiedArray[3*i+2] = 0;

    rectifiedArray[r->N*3+3*i] = 0;
    rectifiedArray[r->N*3+3*i+1] = 0;
    rectifiedArray[r->N*3+3*i+2] = 0;

    dQ_norm = NORM(dQ, i) / NORM(dhem->Qosc, i);

    if(t > dhem->rectifyTimeArray[i] ||
      dQ_norm > r->ri_tes.dq_max)
    {
      rectifyFlag = 1;
      break;
    }
  }

  for(int32_t i = 1; i < r->N; i++)
  {
    if(rectifyFlag != 0)
    {
      rectifiedCount++;

      for(uint32_t j = 0; j < 3; j++)
      {
        double temp_cs = 0;

        Q[3*i+j] = dhem->Qosc[3*i+j];
        add_cs(&Q[3*i+j], &temp_cs, dQ[3*i+j]);
        add_cs(&Q[3*i+j], &temp_cs, r->ri_tes.uVars->uv_csq[3*i+j]);            
        add_cs(&Q[3*i+j], &temp_cs, r->ri_tes.radau->cs_dq[3*i+j]);

        dQ[3*i+j] = -temp_cs;
        r->ri_tes.radau->cs_dq[3*i+j] = 0;
        r->ri_tes.uVars->uv_csq[3*i+j] = 0;
      }

      for(uint32_t j = 0; j < 3; j++)
      {
        double temp_cs = 0;

        P[3*i+j] = dhem->Posc[3*i+j];
        add_cs(&P[3*i+j], &temp_cs, dP[3*i+j]);
        add_cs(&P[3*i+j], &temp_cs, r->ri_tes.uVars->uv_csp[3*i+j]);            
        add_cs(&P[3*i+j], &temp_cs, r->ri_tes.radau->cs_dp[3*i+j]);

        dP[3*i+j] = -temp_cs;
        r->ri_tes.radau->cs_dp[3*i+j] = 0;
        r->ri_tes.uVars->uv_csp[3*i+j] = 0;
      }          

      reb_rebasis_osc_orbits(r, Q, P, t, i);
      dhem->rectifyTimeArray[i] = t + dhem->rectificationPeriod[i];

      rectifiedArray[3*i] = 1;
      rectifiedArray[3*i+1] = 1;
      rectifiedArray[3*i+2] = 1;
      rectifiedArray[r->N*3+3*i] = 1;
      rectifiedArray[r->N*3+3*i+1] = 1;
      rectifiedArray[r->N*3+3*i+2] = 1;
    }    
  }


  return rectifiedCount;
}

static void reb_dhem_calc_osc_orbit_derivs(struct reb_simulation* r, double const * const __restrict__ Qosc, double const * const __restrict__ Posc, 
                                                      double * const __restrict__ Qosc_dot, double * const __restrict__ Posc_dot)
{
  DHEM * dhem = r->ri_tes.rhs;
  const double GM0 = -r->G*dhem->m[0];

  for(int32_t i = 1; i < r->N; i++)
  {
    const double m = r->ri_tes.mass[i];
    const double GMM = GM0*m;
    // Copy across the velocity term.
    Qosc_dot[3*i+0] = Posc[3*i] / m;
    Qosc_dot[3*i+1] = Posc[3*i+1] / m;
    Qosc_dot[3*i+2] = Posc[3*i+2] / m;

    // Calculate the change in momenta term.
    double Q = NORM(Qosc, i);
    double const GMM_Q3 =  GMM/(Q*Q*Q);
    Posc_dot[3*i+0] = GMM_Q3*Qosc[3*i];
    Posc_dot[3*i+1] = GMM_Q3*Qosc[3*i+1];
    Posc_dot[3*i+2] = GMM_Q3*Qosc[3*i+2];
  }
}

static void reb_dhem_calc_osc_orbits_for_all_stages(struct reb_simulation* r, double t0, double h, const double * hArr, uint32_t z_stagesPerStep, uint32_t z_rebasis)
{
  DHEM * dhem = r->ri_tes.rhs;
  if(z_rebasis != 0)
  {
    reb_calc_osc_orbits(r, dhem->XoscArr, t0, h, hArr, z_stagesPerStep, z_rebasis);  
  }
  else
  {
    reb_calc_osc_orbits(r, dhem->XoscPredArr, t0, h, hArr, z_stagesPerStep, z_rebasis);  
  }

  for(uint32_t i = 0; i < z_stagesPerStep; i++)
  {
      double const * const __restrict__ Qout = dhem->XoscArr[i];
      double const * const __restrict__ Pout = &dhem->XoscArr[i][3*r->N];
      double * const __restrict__ Q_dot_out = dhem->Xosc_dotArr[i];
      double * const __restrict__ P_dot_out = &dhem->Xosc_dotArr[i][3*r->N];

      reb_dhem_calc_osc_orbit_derivs(r, Qout, Pout, Q_dot_out, P_dot_out);
  }
}

static double reb_find_min_particle_period(struct reb_simulation* r)
{
    uint32_t N = r->N;  
    double period_min = 1E15;

    for(uint32_t i=1;i<N;i++) 
    {
        double period = r->ri_tes.uVars->period[i];
        period_min = period < period_min ? period : period_min;
    }                        
    
  return period_min;
}

static void reb_dhem_init(struct reb_simulation* r, double z_rectificationPeriodDefault, uint32_t z_stagesPerStep)
{
  // Get memory for the dhem state vectors.
  r->ri_tes.rhs = (DHEM*)malloc(sizeof(DHEM));
  memset(r->ri_tes.rhs, 0, sizeof(DHEM));
  DHEM * dhem = r->ri_tes.rhs;

  dhem->X = (double*)malloc(r->ri_tes.stateVectorSize);
  dhem->rectifyTimeArray = (double*)malloc(r->ri_tes.controlVectorSize);
  dhem->rectificationPeriod = (double*)malloc(r->ri_tes.controlVectorSize);

  // Create space to allow for all of the osculating orbits for a step to be stored.
  dhem->XoscStore = (double*)malloc(z_stagesPerStep*r->ri_tes.stateVectorSize);
  dhem->XoscArr = (double **)malloc(z_stagesPerStep*sizeof(double*));
  dhem->Xosc_dotStore = (double*)malloc(z_stagesPerStep*r->ri_tes.stateVectorSize);
  dhem->Xosc_dotArr = (double **)malloc(z_stagesPerStep*sizeof(double*));

  // Create space to allow for all of the osculating orbits for a step to be stored.
  dhem->XoscPredStore = (double*)malloc(z_stagesPerStep*r->ri_tes.stateVectorSize);
  dhem->XoscPredArr = (double **)malloc(z_stagesPerStep*sizeof(double*));


  // Creat space for osculating orbit compensated summation variables
  dhem->XoscStore_cs = (double*)malloc(z_stagesPerStep*r->ri_tes.stateVectorSize);
  dhem->XoscArr_cs = (double **)malloc(z_stagesPerStep*sizeof(double*));
  
  // Set required arrays to zero.
  memset(dhem->X, 0, r->ri_tes.stateVectorSize);
  memset(dhem->rectifyTimeArray, 0, r->ri_tes.controlVectorSize);
  memset(dhem->rectificationPeriod, 0, r->ri_tes.controlVectorSize);

  memset(dhem->XoscStore, 0, z_stagesPerStep*r->ri_tes.stateVectorSize);
  memset(dhem->XoscArr, 0, z_stagesPerStep*sizeof(double *));
  memset(dhem->Xosc_dotStore, 0, z_stagesPerStep*r->ri_tes.stateVectorSize);
  memset(dhem->Xosc_dotArr, 0, z_stagesPerStep*sizeof(double *));

  memset(dhem->XoscPredStore, 0, z_stagesPerStep*r->ri_tes.stateVectorSize);
  memset(dhem->XoscPredArr, 0, z_stagesPerStep*sizeof(double *));

  memset(dhem->XoscStore_cs, 0, z_stagesPerStep*r->ri_tes.stateVectorSize);
  memset(dhem->XoscArr_cs, 0, z_stagesPerStep*sizeof(double *));

  // To enable easier access to the osculating orbits.
  for(uint32_t i = 0; i < z_stagesPerStep; i++)
  {
    dhem->XoscArr[i] = &dhem->XoscStore[i*r->ri_tes.stateVectorLength];
    dhem->XoscPredArr[i] = &dhem->XoscPredStore[i*r->ri_tes.stateVectorLength];
    dhem->Xosc_dotArr[i] = &dhem->Xosc_dotStore[i*r->ri_tes.stateVectorLength];

    dhem->XoscArr_cs[i] = &dhem->XoscStore_cs[i*r->ri_tes.stateVectorLength];
  }

  // Setup pointers for more human readable access.
  dhem->Qosc = dhem->XoscArr[0];
  dhem->Posc = &dhem->Qosc[3*r->N];

  dhem->Qosc_cs = dhem->XoscArr_cs[0];
  dhem->Posc_cs = &dhem->Qosc_cs[3*r->N];

  dhem->Qosc_dot = dhem->Xosc_dotArr[0];
  dhem->Posc_dot = &dhem->Qosc_dot[3*r->N];

  dhem->Q = dhem->X;
  dhem->P = &dhem->X[3*r->N];

  dhem->m = r->ri_tes.mass;
  dhem->mTotal = 0;

  dhem->m_inv = (double*)malloc(r->N*sizeof(double));

  for(int32_t i = 0; i < r->N; i++)
  {
    dhem->m_inv[i] = 1.0 / dhem->m[i];
    dhem->mTotal += dhem->m[i];
    dhem->rectifyTimeArray[i] = r->t + z_rectificationPeriodDefault;
    dhem->rectificationPeriod[i] = z_rectificationPeriodDefault;
  }
}

static void reb_dhem_free(struct reb_simulation* r)
{
  DHEM * rhs = r->ri_tes.rhs;
  free(rhs->X);
  free(rhs->rectifyTimeArray);
  free(rhs->rectificationPeriod);
  free(rhs->XoscStore);
  free(rhs->XoscArr);
  free(rhs->Xosc_dotStore);
  free(rhs->Xosc_dotArr);
  free(rhs->XoscPredStore);
  free(rhs->XoscPredArr);
  free(rhs->XoscStore_cs);
  free(rhs->XoscArr_cs);
  free(rhs->m_inv);
  free(rhs);
  rhs = NULL;
}

static inline void add_cs(double* out, double* cs, double inp)
{
    const double y = inp - cs[0];
    const double t = out[0] + y;
    cs[0] = (t - out[0]) - y;
    out[0] = t;
}


///////////////////////////////////////////////////////////////////////////////////
// Original TES file: dhem.c
///////////////////////////////////////////////////////////////////////////////////
double reb_single_step(struct reb_simulation* r, double z_t, double dt, double dt_last_done)
{
    RADAU * radau = r->ri_tes.radau;
    uint32_t iterations = 0;
    double dt_new = 0.0;
    uint32_t rectificationCount = reb_dhem_rectify(r, z_t, r->ri_tes.Q_dh, r->ri_tes.P_dh, radau->dQ,
                                        radau->dP, radau->rectifiedArray, FINAL_STAGE_INDEX);
    radau->rectifications += rectificationCount;

    // Calculate the osculating orbits.
    reb_dhem_calc_osc_orbits_for_all_stages(r, z_t, dt, hArr, OSCULATING_ORBIT_SLOTS, 1);

    reb_clear_rectified_b_fields(r, &radau->B, radau->rectifiedArray);
    reb_clear_rectified_b_fields(r, &radau->B_1st, radau->rectifiedArray);

    reb_calc_g_from_b(r); 

    reb_radau_step(r, &iterations, z_t, dt);
    radau->convergenceIterations += iterations;

    dt_new = r->ri_tes.epsilon > 0 ? reb_calc_stepsize(r, dt, dt_last_done, z_t) : dt;

    reb_analytical_continuation(r, &radau->B_1st, &radau->Blast_1st, dt, dt_new, radau->rectifiedArray);
    reb_analytical_continuation(r, &radau->B, &radau->Blast, dt, dt_new, radau->rectifiedArray);

    // Perform a linear update to the drift of the COM.
    for(uint32_t i = 0; i < 3; i++)
    {
      r->ri_tes.COM[i] += r->ri_tes.COM_dot[i]*dt;
    }

    return dt_new;
}

static void reb_radau_init(struct reb_simulation* r)
{
  r->ri_tes.radau = (RADAU *)malloc(sizeof(RADAU));
  memset(r->ri_tes.radau, 0, sizeof(RADAU));
  RADAU * radau = r->ri_tes.radau;

  radau->dX = (double*)malloc(r->ri_tes.stateVectorSize);
  radau->Xout = (double*)malloc(r->ri_tes.stateVectorSize);
  radau->predictors = (double*)malloc(r->ri_tes.stateVectorSize);

  memset(radau->dX, 0, r->ri_tes.stateVectorSize);
  memset(radau->Xout, 0, r->ri_tes.stateVectorSize);
  memset(radau->predictors, 0, r->ri_tes.stateVectorSize);

  radau->dQ = radau->dX;
  radau->dP = &radau->dX[3*r->N];
  radau->Qout = radau->Xout;
  radau->Pout = &radau->Xout[r->ri_tes.stateVectorLength/2];

  // Copy to here so that we are ready to output to a file before we calculate osculating orbtis.
  memcpy(radau->Qout, r->ri_tes.Q_dh, r->ri_tes.stateVectorSize / 2);
  memcpy(radau->Pout, r->ri_tes.P_dh, r->ri_tes.stateVectorSize / 2);

  radau->rectifiedArray = (uint32_t*)malloc(sizeof(uint32_t)*r->ri_tes.stateVectorLength);
  memset(radau->rectifiedArray, 0, sizeof(uint32_t)*r->ri_tes.stateVectorLength);

  radau->b6_store = (double*)malloc(r->ri_tes.stateVectorSize);

  memset(radau->b6_store, 0, r->ri_tes.stateVectorSize);

  //@todo should be able to remove these, but test.
  radau->fCalls = 0;
  radau->rectifications = 0;
  radau->convergenceIterations = 0;

  reb_init_radau_step(r);
}

static void reb_radau_free(struct reb_simulation* r)
{
  RADAU * radau = r->ri_tes.radau;
  reb_free_radau_step(r);
  free(radau->dX);
  free(radau->Xout);
  free(radau->predictors);
  free(radau->rectifiedArray);
  free(radau->b6_store);
  free(radau);
  radau = NULL;
}

double reb_calc_stepsize(struct reb_simulation* r, double h, double hLast, double t)
{
  double hTrial = 0.0;

  // Get the error estimate and orbit size estimate.
  double errMax = reb_calc_step_error(r, h, t);

  if(errMax > r->ri_tes.epsilon*10.0 && r->steps_done == 0)
  {
    printf("\nWarning! The initial step size chosen does not capture the system dynamics. Either do not set the value of sim.dt initially, or if initialising a close encounter then set sim.dt to a smaller initial step size. The current initial step size is %e.", r->dt);
  }
    
  if(isnormal(errMax))
  {
    hTrial = h*pow(r->ri_tes.epsilon / errMax, (1.0/7.0));
  }
  else
  {
    hTrial = 1.1*h;
  }

  // Impose a minimum step size.
  hTrial = hTrial < MIN_STEP_SIZE ? MIN_STEP_SIZE : hTrial;

  // Limit step size growth to 4x per step.
  hTrial = (hTrial > MAX_STEP_SIZE_GROWTH*h) ? h*MAX_STEP_SIZE_GROWTH : hTrial;

  return hTrial;
}


static void reb_clear_rectified_b_fields(struct reb_simulation* r, controlVars * B, uint32_t * rectifiedArray)
{
  for(uint32_t i = 0; i < r->ri_tes.stateVectorLength; i++)
  {
    if(rectifiedArray[i] > 0)
    {
      B->p0[i] = 0.0;
      B->p1[i] = 0.0;
      B->p2[i] = 0.0;
      B->p3[i] = 0.0;
      B->p4[i] = 0.0;
      B->p5[i] = 0.0;
      B->p6[i] = 0.0;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////
// Original TES file: radau_step.c
///////////////////////////////////////////////////////////////////////////////////
static void reb_radau_step(struct reb_simulation* r, uint32_t * z_iterations, double t, double h)
{
  RADAU * radau = r->ri_tes.radau;
  double errMax = 0;
  controlVars * z_csB = &radau->cs_B;
  const uint32_t n3 = 3*r->N;
  controlVars * G = &r->ri_tes.radau->G;
  controlVars * B = &r->ri_tes.radau->B;
  controlVars * G_1st = &r->ri_tes.radau->G_1st;
  controlVars * B_1st = &r->ri_tes.radau->B_1st;

  reb_dhem_rhs_wrapped(r, radau->dQ, radau->dP, radau->dState0, &radau->dState0[3*r->N],
             radau->ddState0, &radau->ddState0[3*r->N], 0, radau->cs_dState0, radau->cs_ddState0);

  radau->fCalls++;

  // Wipe the cs variables for B summation.
  reb_clear_controlvars(z_csB);
  reb_clear_controlvars(&radau->cs_B1st);

  uint32_t n = 1;
  for(n = 1; n < MAX_ITERATIONS; n++)
  {
    z_iterations[0] = n + 1;

    for(uint32_t i = 1; i <= stages; i++)
    {
      reb_calc_predictors(h, hArr[i], radau->dX, radau->dState0, radau->ddState0, B,
                          radau->predictors, radau->cs_dX, 3, r->ri_tes.stateVectorLength/2);
      reb_calc_predictors_1st_order(h, hArr[i], radau->dX, radau->dState0, B_1st, radau->predictors, 
                                  radau->cs_dX, (int)r->ri_tes.stateVectorLength/2, r->ri_tes.stateVectorLength);

      reb_dhem_rhs_wrapped(r, radau->predictors, &radau->predictors[3*r->N], radau->dState, &radau->dState[3*r->N],
           radau->ddState, &radau->ddState[3*r->N], i,
           radau->cs_dState, radau->cs_ddState);

      radau->fCalls++;

      switch(i)
      {
        case 1: 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p0[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];                  
                  G->p0[i] = Gi*rr_inv[0];
                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), G->p0[i]-tmp); 
                }                

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p0[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st->p0[j] = Gi2*rr_inv[0];
                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), G_1st->p0[j]-tmp2);
                }                       
                break;
        case 2: 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p1[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G->p1[i] = (Gi*rr_inv[1] - G->p0[i]) *rr_inv[2];
                  tmp = G->p1[i] - tmp;
                  
                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[0]);
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp);
                }    

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p1[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st->p1[j] = (Gi2*rr_inv[1] - G_1st->p0[j]) *rr_inv[2];
                  tmp2 = G_1st->p1[j] - tmp2;

                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[0]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2);                    
                }                
                break;
        case 3:                 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p2[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  G->p2[i] = ((Gi*rr_inv[3] - G->p0[i]) *rr_inv[4] - G->p1[i])*rr_inv[5];

                  tmp = G->p2[i] - tmp;

                  // We only need to calculate the change in B.
                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[1]);
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp*c[2]);
                  add_cs(&(B->p2[i]), &(radau->cs_B.p2[i]), tmp);            
                }       

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p2[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st->p2[j] = ((Gi2*rr_inv[3] - G_1st->p0[j]) *rr_inv[4] - G_1st->p1[j])*rr_inv[5];

                  tmp2 = G_1st->p2[j] - tmp2;

                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[1]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[2]);
                  add_cs(&(B_1st->p2[j]), &(radau->cs_B1st.p2[j]), tmp2);      
                }                         
                break;
        case 4:  
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p3[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];

                  G->p3[i] = (((Gi*rr_inv[6] - G->p0[i]) *rr_inv[7] - G->p1[i])*rr_inv[8] - G->p2[i]) *rr_inv[9];

                  tmp = G->p3[i] - tmp;

                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[3]);
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp*c[4]);
                  add_cs(&(B->p2[i]), &(radau->cs_B.p2[i]), tmp*c[5]);
                  add_cs(&(B->p3[i]), &(radau->cs_B.p3[i]), tmp);              
                }     

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p3[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];

                  G_1st->p3[j] = (((Gi2*rr_inv[6] - G_1st->p0[j]) *rr_inv[7] - G_1st->p1[j])*rr_inv[8] - G_1st->p2[j]) *rr_inv[9];

                  tmp2 = G_1st->p3[j] - tmp2;

                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[3]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[4]);
                  add_cs(&(B_1st->p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[5]);
                  add_cs(&(B_1st->p3[j]), &(radau->cs_B1st.p3[j]), tmp2);    
                }                

                break;
        case 5:  
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p4[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];                  
                  G->p4[i] = ((((Gi*rr_inv[10] - G->p0[i]) *rr_inv[11] - G->p1[i])*rr_inv[12] - G->p2[i])*rr_inv[13] - G->p3[i]) *rr_inv[14];
                  tmp = G->p4[i] - tmp;
                  

                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[6]);
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp*c[7]);
                  add_cs(&(B->p2[i]), &(radau->cs_B.p2[i]), tmp*c[8]);
                  add_cs(&(B->p3[i]), &(radau->cs_B.p3[i]), tmp*c[9]);
                  add_cs(&(B->p4[i]), &(radau->cs_B.p4[i]), tmp);
                }

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p4[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st->p4[j] = ((((Gi2*rr_inv[10] - G_1st->p0[j]) *rr_inv[11] - G_1st->p1[j])*rr_inv[12] - G_1st->p2[j])*rr_inv[13] - G_1st->p3[j]) *rr_inv[14];
                  tmp2 = G_1st->p4[j] - tmp2;
                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[6]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[7]);
                  add_cs(&(B_1st->p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[8]);
                  add_cs(&(B_1st->p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[9]);
                  add_cs(&(B_1st->p4[j]), &(radau->cs_B1st.p4[j]), tmp2);          
                }                

                break;
        case 6:  
                #pragma GCC ivdep 
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p5[i];

                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G->p5[i] = (((((Gi*rr_inv[15] - G->p0[i]) *rr_inv[16] - G->p1[i])*rr_inv[17] -
                                    G->p2[i]) *rr_inv[18] - G->p3[i]) *rr_inv[19] - G->p4[i]) *rr_inv[20];
                  
                  tmp = G->p5[i] - tmp;
                  

                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[10]);                  
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp*c[11]);
                  add_cs(&(B->p2[i]), &(radau->cs_B.p2[i]), tmp*c[12]);
                  add_cs(&(B->p3[i]), &(radau->cs_B.p3[i]), tmp*c[13]);
                  add_cs(&(B->p4[i]), &(radau->cs_B.p4[i]), tmp*c[14]);
                  add_cs(&(B->p5[i]), &(radau->cs_B.p5[i]), tmp);
                }

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p5[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j]; 
                  G_1st->p5[j] = (((((Gi2*rr_inv[15] - G_1st->p0[j]) *rr_inv[16] - G_1st->p1[j])*rr_inv[17] -
                                    G_1st->p2[j]) *rr_inv[18] - G_1st->p3[j]) *rr_inv[19] - G_1st->p4[j]) *rr_inv[20]; 
                  tmp2 = G_1st->p5[j] - tmp2;

                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[10]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[11]);
                  add_cs(&(B_1st->p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[12]);
                  add_cs(&(B_1st->p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[13]);
                  add_cs(&(B_1st->p4[j]), &(radau->cs_B1st.p4[j]), tmp2*c[14]);
                  add_cs(&(B_1st->p5[j]), &(radau->cs_B1st.p5[j]), tmp2);

                }
                break;
        case 7:  
                radau->acc_ptr = radau->ddState;

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st->p6[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];

                  G_1st->p6[j] = ((((((Gi2*rr_inv[21] - G_1st->p0[j]) *rr_inv[22] - G_1st->p1[j])*rr_inv[23] - G_1st->p2[j]) *rr_inv[24] - 
                                  G_1st->p3[j]) *rr_inv[25] - G_1st->p4[j]) *rr_inv[26] - G_1st->p5[j]) *rr_inv[27];
                  tmp2 = G_1st->p6[j] - tmp2;                                  

                  add_cs(&(B_1st->p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[15]);
                  add_cs(&(B_1st->p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[16]);
                  add_cs(&(B_1st->p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[17]);
                  add_cs(&(B_1st->p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[18]);
                  add_cs(&(B_1st->p4[j]), &(radau->cs_B1st.p4[j]), tmp2*c[19]);
                  add_cs(&(B_1st->p5[j]), &(radau->cs_B1st.p5[j]), tmp2*c[20]);
                  add_cs(&(B_1st->p6[j]), &(radau->cs_B1st.p6[j]), tmp2);
                }

                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G->p6[i]; 
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G->p6[i] =     ((((((Gi*rr_inv[21] - G->p0[i]) *rr_inv[22] - G->p1[i])*rr_inv[23] - G->p2[i]) *rr_inv[24] - 
                                  G->p3[i]) *rr_inv[25] - G->p4[i]) *rr_inv[26] - G->p5[i]) *rr_inv[27];
                                  
                  tmp = G->p6[i] - tmp;
                  
                  add_cs(&(B->p0[i]), &(radau->cs_B.p0[i]), tmp*c[15]);
                  add_cs(&(B->p1[i]), &(radau->cs_B.p1[i]), tmp*c[16]);
                  add_cs(&(B->p2[i]), &(radau->cs_B.p2[i]), tmp*c[17]);
                  add_cs(&(B->p3[i]), &(radau->cs_B.p3[i]), tmp*c[18]);
                  add_cs(&(B->p4[i]), &(radau->cs_B.p4[i]), tmp*c[19]);
                  add_cs(&(B->p5[i]), &(radau->cs_B.p5[i]), tmp*c[20]);
                  add_cs(&(B->p6[i]), &(radau->cs_B.p6[i]), tmp);
                  // Store b6 for convergence criteria.
                  radau->b6_store[i] = tmp;    
                }                
                break;
      }
    }
    double b6Max = 0;
    double accMax = 0;
    double acc_max_tb = 0;
    double ratioMax = 0;
    double ratio = 0;
    errMax = 0;

    // Determine error as compared to two body acceleration
    b6Max = 0;
    accMax = 0;
    acc_max_tb = 0;
    errMax = 0;
    double q_ddot[3*r->N];

    for(int32_t j = 1; j < r->N; j++)
    {
      q_ddot[3*j+0] = r->ri_tes.rhs->Xosc_dotArr[8][3*r->N+3*j+0] / r->ri_tes.mass[j];
      q_ddot[3*j+1] = r->ri_tes.rhs->Xosc_dotArr[8][3*r->N+3*j+1] / r->ri_tes.mass[j];
      q_ddot[3*j+2] = r->ri_tes.rhs->Xosc_dotArr[8][3*r->N+3*j+2] / r->ri_tes.mass[j];
    }    

    for(int32_t j = 3; j < 3*r->N; j++)
    {
      double b6 = fabs(radau->b6_store[j]);
      double acc = fabs(radau->acc_ptr[j]);

      // double acc_tb = fabs(sim->rhs->Xosc_dotArr[8][j]);
      double acc_tb = fabs(q_ddot[j]);

      if(isnormal(acc_tb))
      {
        acc_max_tb = acc_tb > acc_max_tb ? acc_tb : acc_max_tb;
      }

      if(isnormal(b6))
      {
        b6Max = b6 > b6Max ? b6 : b6Max;
      }

      if(isnormal(acc))
      {
        accMax = acc > accMax ? acc : accMax;
      }
    }
    errMax = b6Max/acc_max_tb;
  

    for(int32_t k = 1; k < r->N; k++)
    {
      ratio = NORM(radau->dX, k) / NORM(r->ri_tes.rhs->Qosc, k);

      ratioMax = ratio > ratioMax ? ratio : ratioMax;
    }

    // Have we converged?
    if(((errMax < 1E-15 ) && n >= MIN_ITERATIONS))
    {
      z_iterations[0] = n;
      break;
    }
  }
 
  if(n >= MAX_ITERATIONS && MAX_ITERATIONS != MIN_ITERATIONS)
  {
    z_iterations[0] = MAX_ITERATIONS;
  }

  // Apply correctors to obtain X at t=1.
  reb_update_state(h, radau->dState0, radau->ddState0, B, radau->dX, 
                    radau->cs_dX, 0, (int)(n3));

  reb_update_state_1st_order(h, radau->dState0, B_1st, radau->dX, 
                                radau->cs_dX, (int)(n3), r->ri_tes.stateVectorLength);                                

  // Now that we have cs_dX we can perform our correction.
  reb_apply_osc_orbit_corrector(r, r->ri_tes.rhs->XoscArr, t+h, 9); 
}

double reb_calc_step_error(struct reb_simulation* r, double h, double t)
{
  RADAU * radau = r->ri_tes.radau;
  double b6Max = 0;
  double accMax = 0;
  double errMax = 0;

  for(int32_t i = 0; i < r->N; i++) 
  {
      for(uint32_t j = 0; j < 3; j++)
      {
        double b6 = fabs(radau->B.p6[3*i+j]);
        double acc = fabs(radau->acc_ptr[3*i+j]); 

        if(isnormal(b6))
        {
          b6Max = b6 > b6Max ? b6 : b6Max;
        }

        if(isnormal(acc))
        {
          accMax = acc > accMax ? acc : accMax;
        }
      }
  }  
  errMax = b6Max/accMax;

  return errMax;
}

static void reb_calc_g_from_b_internal(controlVars * z_G, controlVars * z_B, uint32_t z_start, uint32_t z_end)
{
  for(uint32_t i = z_start; i < z_end; i++)
  {
      z_G->p0[i] = z_B->p6[i]*d[15] + z_B->p5[i]*d[10] + z_B->p4[i]*d[6] + z_B->p3[i]*d[3]  + z_B->p2[i]*d[1]  + z_B->p1[i]*d[0]  + z_B->p0[i];
      z_G->p1[i] = z_B->p6[i]*d[16] + z_B->p5[i]*d[11] + z_B->p4[i]*d[7] + z_B->p3[i]*d[4]  + z_B->p2[i]*d[2]  + z_B->p1[i];
      z_G->p2[i] = z_B->p6[i]*d[17] + z_B->p5[i]*d[12] + z_B->p4[i]*d[8] + z_B->p3[i]*d[5]  + z_B->p2[i];
      z_G->p3[i] = z_B->p6[i]*d[18] + z_B->p5[i]*d[13] + z_B->p4[i]*d[9] + z_B->p3[i];
      z_G->p4[i] = z_B->p6[i]*d[19] + z_B->p5[i]*d[14] + z_B->p4[i];
      z_G->p5[i] = z_B->p6[i]*d[20] + z_B->p5[i];
      z_G->p6[i] = z_B->p6[i];
  }
}

static void reb_calc_g_from_b(struct reb_simulation* r)
{
  reb_calc_g_from_b_internal(&r->ri_tes.radau->G, &r->ri_tes.radau->B, 0, (int)(r->ri_tes.stateVectorLength/2));
  reb_calc_g_from_b_internal(&r->ri_tes.radau->G_1st, &r->ri_tes.radau->B_1st, (int)(r->ri_tes.stateVectorLength/2), r->ri_tes.stateVectorLength);
}

static void reb_analytical_continuation(struct reb_simulation* r, controlVars * z_B, controlVars * z_Blast, const double h,
                             const double h_new, const uint32_t * const rectificationArray)
{
  const double ratio = h_new / h;
  const double q1 = ratio;
  const double q2 = q1 * q1;
  const double q3 = q1 * q2;
  const double q4 = q2 * q2;
  const double q5 = q2 * q3;
  const double q6 = q3 * q3;
  const double q7 = q3 * q4;

  for(uint32_t i = 0; i < r->ri_tes.stateVectorLength; i++)
  {
    double dB0 = 0;
    double dB1 = 0;
    double dB2 = 0;
    double dB3 = 0;
    double dB4 = 0;
    double dB5 = 0;
    double dB6 = 0;

    // Removed these when trying to resolve the iterations problem after rectification.
    // // If we have rectified then we cant use the delta from last step to improve B.
    // if(rectificationArray == NULL || (rectificationArray != NULL && rectificationArray[i] == 0))
    // {
    //   dB0 = z_B->p0[i] - z_Blast->p0[i];
    //   dB1 = z_B->p1[i] - z_Blast->p1[i];
    //   dB2 = z_B->p2[i] - z_Blast->p2[i];
    //   dB3 = z_B->p3[i] - z_Blast->p3[i];
    //   dB4 = z_B->p4[i] - z_Blast->p4[i];
    //   dB5 = z_B->p5[i] - z_Blast->p5[i];
    //   dB6 = z_B->p6[i] - z_Blast->p6[i];
    // }

    // Save the predicted values of B.
    z_Blast->p0[i] = q1*(z_B->p6[i]* 7.0 + z_B->p5[i]* 6.0 + z_B->p4[i]* 5.0 + z_B->p3[i]* 4.0 + z_B->p2[i]* 3.0 + z_B->p1[i]*2.0 + z_B->p0[i]);
    z_Blast->p1[i] = q2*(z_B->p6[i]*21.0 + z_B->p5[i]*15.0 + z_B->p4[i]*10.0 + z_B->p3[i]* 6.0 + z_B->p2[i]* 3.0 + z_B->p1[i]);
    z_Blast->p2[i] = q3*(z_B->p6[i]*35.0 + z_B->p5[i]*20.0 + z_B->p4[i]*10.0 + z_B->p3[i]* 4.0 + z_B->p2[i]);
    z_Blast->p3[i] = q4*(z_B->p6[i]*35.0 + z_B->p5[i]*15.0 + z_B->p4[i]* 5.0 + z_B->p3[i]);
    z_Blast->p4[i] = q5*(z_B->p6[i]*21.0 + z_B->p5[i]* 6.0 + z_B->p4[i]);
    z_Blast->p5[i] = q6*(z_B->p6[i]* 7.0 + z_B->p5[i]);
    z_Blast->p6[i] = q7* z_B->p6[i];

    if(rectificationArray == NULL || (rectificationArray != NULL && rectificationArray[i] == 0))
    {
      z_B->p0[i] = z_Blast->p0[i] + dB0;
      z_B->p1[i] = z_Blast->p1[i] + dB1;
      z_B->p2[i] = z_Blast->p2[i] + dB2;
      z_B->p3[i] = z_Blast->p3[i] + dB3;
      z_B->p4[i] = z_Blast->p4[i] + dB4;
      z_B->p5[i] = z_Blast->p5[i] + dB5;
      z_B->p6[i] = z_Blast->p6[i] + dB6;
    }    
  }
}

static void reb_calc_predictors_1st_order(double h, double hSample, double * z_state0, double * z_dState0, 
                                    controlVars * z_B, double * z_predictors, double * z_csState, 
                                    uint32_t z_start, uint32_t z_end)
{
    double s[9];
    s[0] = h * hSample;
    s[1] = s[0] * hSample / 2.;
    s[2] = 2. * s[1] * hSample / 3.;
    s[3] = 3. * s[2] * hSample / 4.;
    s[4] = 4. * s[3] * hSample / 5.;
    s[5] = 5. * s[4] * hSample / 6.;
    s[6] = 6. * s[5] * hSample / 7.;
    s[7] = 7. * s[6] * hSample / 8.;

    controlVars_const Bconst = control_vars_cast(z_B[0]);


    #pragma GCC ivdep
    for(uint32_t i = z_start; i < z_end; i++)
    {
      double pred = (s[7]*Bconst.p6[i] + s[6]*Bconst.p5[i] + s[5]*Bconst.p4[i] + 
                    s[4]*Bconst.p3[i] + s[3]*Bconst.p2[i] + s[2]*Bconst.p1[i] + s[1]*Bconst.p0[i] + 
                    s[0]*z_dState0[i]);
      z_predictors[i] = pred + z_state0[i];            
    }
}

            
static void reb_calc_predictors(double h, double hSample, double const * __restrict__ z_state0, double const * __restrict__ z_dState, 
                         double const * __restrict__ z_ddState, controlVars const * z_B, double * __restrict__ z_predictors, 
                         double const * __restrict__ z_csState, uint32_t const z_start, uint32_t const z_end)
{
    double s[9];
    s[0] = h * hSample;
    s[1] = s[0] * s[0] / 2.;
    s[2] = s[1] * hSample / 3.;
    s[3] = s[2] * hSample / 2.;
    s[4] = 3. * s[3] * hSample / 5.;
    s[5] = 2. * s[4] * hSample / 3.;
    s[6] = 5. * s[5] * hSample / 7.;
    s[7] = 3. * s[6] * hSample / 4.;
    s[8] = 7. * s[7] * hSample / 9.;

    controlVars_const Bconst = control_vars_cast(z_B[0]);

    #pragma GCC ivdep
    for(uint32_t i = z_start; i < z_end; i++)
    {
      double pred = (s[8]*Bconst.p6[i] + s[7]*Bconst.p5[i] + s[6]*Bconst.p4[i] + 
                    s[5]*Bconst.p3[i] + s[4]*Bconst.p2[i] + s[3]*Bconst.p1[i] + s[2]*Bconst.p0[i] + 
                    s[1]*z_ddState[i] + s[0]*z_dState[i]);
      z_predictors[i] = pred + z_state0[i];            
    }
}

static void reb_update_state_1st_order(double h, double * z_dState, controlVars * z_B, double * z_state0, 
                                double * z_csState, uint32_t z_start, uint32_t z_end)
{
  for(uint32_t i = z_start; i < z_end; i++)
  {
    add_cs(&z_state0[i], &z_csState[i], z_B->p6[i] / 8.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p5[i] / 7.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p4[i] / 6.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p3[i] / 5.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p2[i] / 4.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p1[i] / 3.0*h);
    add_cs(&z_state0[i], &z_csState[i], z_B->p0[i] / 2.0*h);
    add_cs(&z_state0[i], &z_csState[i], h*z_dState[i]);    
  }
}

static void reb_update_state(double h, double * z_dState, double * z_ddState,
                        controlVars * z_B, double * z_state0, double * z_csState, uint32_t z_start, uint32_t z_end)
{
  double h2 = h*h;
  for(uint32_t i = z_start; i < z_end; i++)
  {
      add_cs(&z_state0[i], &z_csState[i], z_B->p6[i] / 72.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p5[i] / 56.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p4[i] / 42.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p3[i] / 30.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p2[i] / 20.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p1[i] / 12.0 * h2);
      add_cs(&z_state0[i], &z_csState[i], z_B->p0[i] / 6.0 *  h2);
      add_cs(&z_state0[i], &z_csState[i], z_ddState[i] / 2.0 *h2);
      add_cs(&z_state0[i], &z_csState[i], h*z_dState[i]);
  }
}

static void reb_init_radau_step(struct reb_simulation* r)
{
    RADAU * radau = r->ri_tes.radau;
    reb_init_controlvars(&radau->G, r->ri_tes.stateVectorSize);    
    reb_init_controlvars(&radau->B, r->ri_tes.stateVectorSize);    
    reb_init_controlvars(&radau->Blast, r->ri_tes.stateVectorSize);   
    reb_init_controlvars(&radau->Blast_1st, r->ri_tes.stateVectorSize);
    reb_init_controlvars(&radau->G_1st, r->ri_tes.stateVectorSize);
    reb_init_controlvars(&radau->B_1st, r->ri_tes.stateVectorSize);
    
    reb_init_controlvars(&radau->cs_B1st, r->ri_tes.stateVectorSize);
    reb_init_controlvars(&radau->cs_B, r->ri_tes.stateVectorSize);

    radau->dState0 = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->ddState0 = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->dState = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->ddState = (double *)malloc(r->ri_tes.stateVectorSize);

    memset(radau->dState0, 0, r->ri_tes.stateVectorSize);
    memset(radau->ddState0, 0, r->ri_tes.stateVectorSize);
    memset(radau->dState, 0, r->ri_tes.stateVectorSize);
    memset(radau->ddState, 0, r->ri_tes.stateVectorSize);

    // Compensated summation arrays
    radau->cs_dState0 = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->cs_ddState0 = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->cs_dState = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->cs_ddState = (double *)malloc(r->ri_tes.stateVectorSize);

    memset(radau->cs_dState0, 0, r->ri_tes.stateVectorSize);
    memset(radau->cs_ddState0, 0, r->ri_tes.stateVectorSize);
    memset(radau->cs_dState, 0, r->ri_tes.stateVectorSize);
    memset(radau->cs_ddState, 0, r->ri_tes.stateVectorSize);

    radau->cs_dX = (double *)malloc(r->ri_tes.stateVectorSize);
    radau->cs_dq = radau->cs_dX;
    radau->cs_dp = &radau->cs_dX[(int)r->ri_tes.stateVectorLength/2];

    memset(radau->cs_dX, 0, r->ri_tes.stateVectorSize);
}

static void reb_free_radau_step(struct reb_simulation* r)
{
  RADAU * radau = r->ri_tes.radau;
  reb_free_controlvars(&r->ri_tes.radau->G);
  reb_free_controlvars(&r->ri_tes.radau->B);
  reb_free_controlvars(&r->ri_tes.radau->Blast);
  reb_free_controlvars(&r->ri_tes.radau->Blast_1st);
  reb_free_controlvars(&r->ri_tes.radau->G_1st);
  reb_free_controlvars(&r->ri_tes.radau->B_1st);    
  reb_free_controlvars(&r->ri_tes.radau->cs_B1st);
  reb_free_controlvars(&r->ri_tes.radau->cs_B);
  free(radau->dState0);
  free(radau->ddState0);
  free(radau->dState);
  free(radau->ddState);
  free(radau->cs_dState0);
  free(radau->cs_ddState0);
  free(radau->cs_dState);
  free(radau->cs_ddState);
  free(radau->cs_dX);
}


static void reb_init_controlvars(controlVars * var, uint32_t size)
{
  var->size = size;
  var->p0 = (double*)malloc(size);
  var->p1 = (double*)malloc(size);
  var->p2 = (double*)malloc(size);
  var->p3 = (double*)malloc(size);
  var->p4 = (double*)malloc(size);
  var->p5 = (double*)malloc(size);
  var->p6 = (double*)malloc(size);

  memset(var->p0, 0, size);
  memset(var->p1, 0, size);
  memset(var->p2, 0, size);
  memset(var->p3, 0, size);
  memset(var->p4, 0, size);
  memset(var->p5, 0, size);
  memset(var->p6, 0, size);
}

static void reb_free_controlvars(controlVars * var)
{
  free(var->p0);
  free(var->p1);
  free(var->p2);
  free(var->p3);
  free(var->p4);
  free(var->p5);
  free(var->p6);
}

static void reb_clear_controlvars(controlVars * var)
{
  for(uint32_t i = 0; i < var->size/sizeof(double); i++)
  {
    var->p0[i] = 0.0;
    var->p1[i] = 0.0;
    var->p2[i] = 0.0;
    var->p3[i] = 0.0;
    var->p4[i] = 0.0;
    var->p5[i] = 0.0;
    var->p6[i] = 0.0;
  }
}

///////////////////////////////////////////////////////////////////////////////////
// Original TES file: UniversalVars.c
///////////////////////////////////////////////////////////////////////////////////
static void reb_calc_osc_orbits(struct reb_simulation* r, double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis)
{
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
  double t_last_rebasis = 0;

  for(uint32_t stage = 0; stage < z_stagesPerStep; stage++)
  {
    double t = t0 + h*h_array[stage];
    
    double dt = 0;

    // This will not allow us to maintain the osculating orbit for more than 
    // a step as we would need a dt = t_rebasis + h*(h_arry[] ...) version.
    if(stage == 0)
    {
      dt = 0;
    }
    else
    {
      dt = h*(h_array[stage]-t_last_rebasis); 
    }

    double C[4] = {0.0, 0.0, 0.0, 0.0};
    for(int32_t i = 1; i < r->N; i++)
    {  
      // Calculate the dt value and wrap around the orbital period.
      dt = fmod(dt, p_uVars->period[i]);

      // Calculate our step since last time we were called and update storage of tLast.
      double h = t - p_uVars->tLast[i];
      p_uVars->tLast[i] = t;

      reb_solve_for_universal_anomaly(r, dt, h, i, C);

      p_uVars->C.c0[i] = C[0];
      p_uVars->C.c1[i] = C[1];
      p_uVars->C.c2[i] = C[2];
      p_uVars->C.c3[i] = C[3];
    }

    for(int32_t i = 1; i < r->N; i++)
    {    
      double * Qout = Xosc_map[stage];
      double * Pout = &Qout[3*r->N];      
      const double X = p_uVars->X[i];
      const double X2 = X*X;
      const double X3 = X2*X;
      const double g1 = p_uVars->C.c1[i]*X;
      const double g2 = p_uVars->C.c2[i]*X2;
      const double g3 = p_uVars->C.c3[i]*X3;

      const double mu = p_uVars->mu;

      // Calculate G
      double dt_h = dt;
      double dt_t = dt-dt_h;
      add_cs(&dt_h, &dt_t, -mu*g3);
      double g = dt_h;
      const double f = -p_uVars->mu*g2 / p_uVars->Q0_norm[i];    

      double Q0x_h = p_uVars->Q0[i*3];
      double Q0x_t = p_uVars->uv_csq[3*i];
      double Q0y_h = p_uVars->Q0[i*3+1];
      double Q0y_t = p_uVars->uv_csq[3*i+1];     
      double Q0z_h = p_uVars->Q0[i*3+2];
      double Q0z_t = p_uVars->uv_csq[3*i+2];;

      double V0x_h = p_uVars->V0[i*3];
      double V0x_t = p_uVars->uv_csv[3*i];
      double V0y_h = p_uVars->V0[i*3+1];
      double V0y_t = p_uVars->uv_csv[3*i+1];
      double V0z_h = p_uVars->V0[i*3+2];
      double V0z_t = p_uVars->uv_csv[3*i+2];

      double dx, dy, dz = 0;
      dx = (((g*V0x_t+f*Q0x_t) + g*V0x_h) + f*Q0x_h);
      dy = (((g*V0y_t+f*Q0y_t) + g*V0y_h) + f*Q0y_h);
      dz = (((g*V0z_t+f*Q0z_t) + g*V0z_h) + f*Q0z_h);

      // Set up summation vars
      Qout[i*3]   = p_uVars->Q0[i*3];
      Qout[i*3+1] = p_uVars->Q0[i*3+1];
      Qout[i*3+2] = p_uVars->Q0[i*3+2];

      // Update from f+g functions.
      add_cs(&Qout[i*3]  , &p_uVars->uv_csq[3*i], dx);
      add_cs(&Qout[i*3+1], &p_uVars->uv_csq[3*i+1], dy);
      add_cs(&Qout[i*3+2], &p_uVars->uv_csq[3*i+2], dz);

      // Copy everything back into storage
      p_uVars->Q1[3*i]   = Qout[i*3];
      p_uVars->Q1[3*i+1] = Qout[i*3+1];
      p_uVars->Q1[3*i+2] = Qout[i*3+2];
      
      const double Qnorm = NORM(p_uVars->Q1, i);
      const double fp = -p_uVars->mu * g1 / (p_uVars->Q0_norm[i]*Qnorm);
      const double gp =  -p_uVars->mu * g2 / Qnorm;

      const double mi =  r->ri_tes.mass[i]; 
      dx = (((gp*V0x_t+fp*Q0x_t) + gp*V0x_h) + fp*Q0x_h)*mi;
      dy = (((gp*V0y_t+fp*Q0y_t) + gp*V0y_h) + fp*Q0y_h)*mi;
      dz = (((gp*V0z_t+fp*Q0z_t) + gp*V0z_h) + fp*Q0z_h)*mi;


      Pout[i*3]   = p_uVars->P0[i*3];
      Pout[i*3+1] = p_uVars->P0[i*3+1];
      Pout[i*3+2] = p_uVars->P0[i*3+2];

      add_cs(&Pout[3*i], &p_uVars->uv_csp[3*i], dx);
      add_cs(&Pout[3*i+1], &p_uVars->uv_csp[3*i+1], dy);
      add_cs(&Pout[3*i+2], &p_uVars->uv_csp[3*i+2], dz);

      // Copy everything back into storage
      p_uVars->P1[3*i]   = Pout[i*3];
      p_uVars->P1[3*i+1] = Pout[i*3+1];
      p_uVars->P1[3*i+2] = Pout[i*3+2];

      if(stage < z_stagesPerStep-1)
      {       
        if(z_rebasis)
        {
          reb_rebasis_osc_orbits(r, Qout, Pout, t, i);
        }
      }
    }
    
    if(z_rebasis != 0)
    {
      t_last_rebasis = h_array[stage];
    }
  } 
}

static void reb_apply_osc_orbit_corrector(struct reb_simulation* r, double **Xosc_map, double t, uint32_t z_stagePerStep)
{ 
    UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
    double * Qout = Xosc_map[z_stagePerStep-1];
    double * Pout = &Qout[3*r->N];      

    for(int32_t i = 1; i < r->N; i++)
    {
        for(uint32_t j = 0; j < 3; j++)
        {
          double csq = 0;
          double csp = 0;

          // Correct the osculating orbit inline with the various cs vars.
          add_cs(&Qout[i*3+j], &csq, -r->ri_tes.radau->cs_dq[3*i+j]);
          add_cs(&Qout[i*3+j], &csq, -p_uVars->uv_csq[3*i+j]);
          add_cs(&Pout[i*3+j], &csp, -r->ri_tes.radau->cs_dp[3*i+j]);
          add_cs(&Pout[i*3+j], &csp, -p_uVars->uv_csp[3*i+j]);   

          p_uVars->uv_csq[3*i+j] = 0;
          p_uVars->uv_csp[3*i+j] = 0;  
          r->ri_tes.radau->cs_dq[3*i+j] = 0;
          r->ri_tes.radau->cs_dp[3*i+j] = 0;        
          
          // Take anything left over from the oscualting orbit corrector and place into the delta.
          add_cs(&r->ri_tes.radau->dQ[3*i+j], &p_uVars->uv_csq[3*i+j], -csq);
          add_cs(&r->ri_tes.radau->dP[3*i+j], &p_uVars->uv_csp[3*i+j], -csp);

          // Obtain cs sum for v for next step.
          p_uVars->uv_csv[3*i+j] = p_uVars->uv_csp[3*i+j]/r->ri_tes.mass[i];
        }
        reb_rebasis_osc_orbits(r, Qout, Pout, t, i);  
    }
}

static void reb_rebasis_osc_orbits(struct reb_simulation* r, double * z_Q, double * z_P, double z_t, uint32_t i)
{
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
  double * mass = r->ri_tes.mass;
  p_uVars->Q0[3*i+0] = z_Q[3*i+0];
  p_uVars->Q0[3*i+1] = z_Q[3*i+1];
  p_uVars->Q0[3*i+2] = z_Q[3*i+2];
  
  p_uVars->V0[3*i+0] = z_P[3*i+0] / mass[i];
  p_uVars->V0[3*i+1] = z_P[3*i+1] / mass[i];
  p_uVars->V0[3*i+2] = z_P[3*i+2] / mass[i];

  p_uVars->uv_csv[3*i] = p_uVars->uv_csp[3*i] / mass[i];
  p_uVars->uv_csv[3*i+1] = p_uVars->uv_csp[3*i+1] / mass[i];
  p_uVars->uv_csv[3*i+2] = p_uVars->uv_csp[3*i+2] / mass[i];

  p_uVars->P0[3*i+0] = z_P[3*i+0];
  p_uVars->P0[3*i+1] = z_P[3*i+1];
  p_uVars->P0[3*i+2] = z_P[3*i+2];

  const double mu = p_uVars->mu;
  const double PIx4_MU = PI_SQ_X4 / mu;

  const double Q0x = z_Q[3*i];
  const double Q0y = z_Q[3*i+1];
  const double Q0z = z_Q[3*i+2];
  const double V0x = z_P[3*i]/mass[i];
  const double V0y = z_P[3*i+1]/mass[i];
  const double V0z = z_P[3*i+2]/mass[i];

  p_uVars->t0[i] = z_t;
  const double Q0_norm = sqrt(Q0x*Q0x+Q0y*Q0y+Q0z*Q0z);
  const double V0_norm = sqrt(V0x*V0x+V0y*V0y+V0z*V0z);
  const double V0_norm2 = V0_norm * V0_norm;
  p_uVars->Q0_norm[i] = Q0_norm;
  
  p_uVars->eta[i] = (double)(Q0x*V0x + Q0y*V0y + Q0z*V0z);
  const double beta = (2.0*mu/Q0_norm)-V0_norm2;
  p_uVars->beta[i] = beta;
  p_uVars->zeta[i] = mu-beta*Q0_norm;

  double a = mu / beta;
  p_uVars->period[i] = sqrt(PIx4_MU*(a*a*a));
  p_uVars->Xperiod[i] = 2.0*PI / sqrt(beta);
  p_uVars->X[i] = 0.0;
  p_uVars->dt[i] = 0;
}


static double reb_solve_for_universal_anomaly(struct reb_simulation* const r, double dt, double h, uint32_t i, double * C)
{
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
  double X = p_uVars->X[i] + (h / p_uVars->Q0_norm[i]);
  X *= (1.0 - X*p_uVars->eta[i]*0.5/p_uVars->Q0_norm[i]);

  X = fmod(X, p_uVars->Xperiod[i]);

  double prevVals[MAX_NEWTON_ITERATIONS];

  for(uint32_t k = 0; k < MAX_NEWTON_ITERATIONS; k++)
  {    
    X = reb_calc_osc_update_value(r, X, dt, i, C);

    for(uint32_t j = 0; j < k; j++)
    {
      if(X == prevVals[j])
      {
        p_uVars->X[i] = X;
        return X;
      }
    }
    prevVals[k] = X;
  }
  return -1.0;
}

static double reb_calc_osc_update_value(struct reb_simulation* const r, double X, double dt, uint32_t i, double * C)
{
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
  const double X2 = X*X;
  const double X3 = X2*X;
  double Z = p_uVars->beta[i]*X2;

  reb_c_stumpff(C, Z);

  const double num = (X*(p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2]) -
              p_uVars->eta[i]*X2*C[2] -
              p_uVars->zeta[i]*X3*C[3] + dt);
  const double denom = p_uVars->Q0_norm[i] + p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2];

  return (num / denom);
}

static void reb_c_stumpff(double * cs, double z)
{
    unsigned int n = 0;
    while(fabs(z) > 0.1){
        z = z/4.0;
        n++;
    }
    
    double c_odd  = invfactorial[STUMPF_ITERATIONS];
    double c_even = invfactorial[STUMPF_ITERATIONS-1];
    for(int np=STUMPF_ITERATIONS-2;np>=3;np-=2){
        c_odd  = invfactorial[np]    - z *c_odd;
        c_even = invfactorial[np-1]  - z *c_even;
    }

    cs[3] = c_odd;
    cs[2] = c_even;
    cs[1] = invfactorial[1]  - z *c_odd;
    cs[0] = invfactorial[0]  - z *c_even;

    for (;n>0;n--){
        cs[3] = (cs[2]+cs[0]*cs[3])/4;
        cs[2] = cs[1]*cs[1]/2;
        cs[1] = cs[0]*cs[1];
        cs[0] = 2.0*cs[0]*cs[0]-1.0;
    }
}

static void reb_init_uvars(struct reb_simulation* const r)
{
  int N = r->N;
  // Create the main control data structure for universal variables.
  r->ri_tes.uVars = (UNIVERSAL_VARS *)malloc(sizeof(UNIVERSAL_VARS));
  memset(r->ri_tes.uVars, 0, sizeof(UNIVERSAL_VARS));
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;

  p_uVars->stateVectorSize = 3 * N * sizeof(double);
  p_uVars->controlVectorSize = N * sizeof(double);

  // Allocate memory for all objects used in universal vars
  p_uVars->Q0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->V0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->Q1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->V1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->P0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->P1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->t0 = (double *)malloc(N*sizeof(double));
  p_uVars->tLast = (double *)malloc(N*sizeof(double));
  p_uVars->Q0_norm = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->beta = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->eta =  (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->zeta = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->period = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->Xperiod = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->X = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->dt = (double *)malloc(p_uVars->controlVectorSize);

  p_uVars->mu = (double)r->G*(double)r->ri_tes.mass[0];
  p_uVars->e = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->a = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->h = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->h_norm = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->peri = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->apo = (double *)malloc(p_uVars->controlVectorSize);
  
  p_uVars->C.c0 = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->C.c1 = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->C.c2 = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->C.c3 = (double *)malloc(p_uVars->controlVectorSize);

  memset(p_uVars->Q0, 0, p_uVars->stateVectorSize);
  memset(p_uVars->V0, 0, p_uVars->stateVectorSize);
  memset(p_uVars->Q1, 0, p_uVars->stateVectorSize);
  memset(p_uVars->V1, 0, p_uVars->stateVectorSize);
  memset(p_uVars->P0, 0, p_uVars->stateVectorSize);
  memset(p_uVars->P1, 0, p_uVars->stateVectorSize);
  memset(p_uVars->t0, 0, N*sizeof(double));
  memset(p_uVars->tLast, 0, N*sizeof(double));
  memset(p_uVars->Q0_norm, 0, p_uVars->controlVectorSize);
  memset(p_uVars->beta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->eta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->zeta, 0, p_uVars->controlVectorSize);
  memset(p_uVars->period, 0, p_uVars->controlVectorSize);
  memset(p_uVars->Xperiod, 0, p_uVars->controlVectorSize);
  memset(p_uVars->X, 0, p_uVars->controlVectorSize);
  memset(p_uVars->a, 0, p_uVars->controlVectorSize);
  memset(p_uVars->e, 0, p_uVars->controlVectorSize);
  memset(p_uVars->h, 0, p_uVars->stateVectorSize);
  memset(p_uVars->h_norm, 0, p_uVars->controlVectorSize);
  memset(p_uVars->peri, 0, p_uVars->controlVectorSize);
  memset(p_uVars->apo, 0, p_uVars->controlVectorSize);
  memset(p_uVars->dt, 0, p_uVars->controlVectorSize);

  // CS vars
  p_uVars->uv_csq = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->uv_csv = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->uv_csp = (double *)malloc(p_uVars->stateVectorSize);
  memset(p_uVars->uv_csq, 0, p_uVars->stateVectorSize);
  memset(p_uVars->uv_csv, 0, p_uVars->stateVectorSize);
  memset(p_uVars->uv_csp, 0, p_uVars->stateVectorSize);  
}

static void reb_free_uvars(struct reb_simulation* r)
{
  UNIVERSAL_VARS * p_uVars = r->ri_tes.uVars;
  free(p_uVars->Q0);
  free(p_uVars->V0);
  free(p_uVars->Q1);
  free(p_uVars->V1);
  free(p_uVars->P0);
  free(p_uVars->P1);
  free(p_uVars->t0);
  free(p_uVars->tLast);
  free(p_uVars->Q0_norm);
  free(p_uVars->beta);
  free(p_uVars->eta);
  free(p_uVars->zeta);
  free(p_uVars->period);
  free(p_uVars->Xperiod);
  free(p_uVars->X);
  free(p_uVars->dt);
  free(p_uVars->e);
  free(p_uVars->a);
  free(p_uVars->h);
  free(p_uVars->h_norm);
  free(p_uVars->peri);
  free(p_uVars->apo);
  free(p_uVars->C.c0);
  free(p_uVars->C.c1);
  free(p_uVars->C.c2);
  free(p_uVars->C.c3);
  free(p_uVars->uv_csq);
  free(p_uVars->uv_csv);
  free(p_uVars->uv_csp);
  free(p_uVars);
}
