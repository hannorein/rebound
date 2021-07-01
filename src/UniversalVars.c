// TES is an open source integration package for modelling exoplanet evolution.
// Copyright (C) <2021>  <Peter Bartram, Alexander Wittig>

// TES is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>. 

#include "UniversalVars.h"
#include "radau.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "dhem.h"

#define PI 3.141592653589793238462643383279
#define PI_SQ_X4 (4.0*PI*PI)
#define PI_X2 (2.0*PI)


double invfactorial[35] = {1.0,1.0,0.5,0.166666666666666666666666666666667,0.0416666666666666666666666666666667,
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

#define NORM(x, i) sqrt((x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2]))
#define DOT(x, y, i) x[3*i+0]*y[3*i+0]+x[3*i+1]*y[3*i+1]+x[3*i+2]*y[3*i+2]

#define CROSS(a,b,c,i) \
	(a)[3*i+0] = (b)[3*i+1] * (c)[3*i+2] - (c)[3*i+1] * (b)[3*i+2]; \
	(a)[3*i+1] = (b)[3*i+2] * (c)[3*i+0] - (c)[3*i+2] * (b)[3*i+0]; \
	(a)[3*i+2] = (b)[3*i+0] * (c)[3*i+1] - (c)[3*i+0] * (b)[3*i+1];

#define CROSS_LOCAL(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define MAX_NEWTON_ITERATIONS 50
#define STUMPF_ITERATIONS 13 // must be an odd number

static UNIVERSAL_VARS * p_uVars = NULL;
static SIMULATION * sim = NULL;

// Internal functions.
static void C_Stumpff(double * cs_in, double z_in);
static double SolveForUnivsersalAnomaly(double dt, double h, uint32_t i, double * C);
static double CalcUpdateValue(double X, double dt, uint32_t i, double * C);
static void CalculateClassicalOrbitalElementsSingle(uint32_t i);
static inline void add_cs(double * out, double * cs, double inp);


void CalculateOsculatingOrbitsForSingleStep(double **Xosc_map, 
                                            const double t0, const double h, double const * const h_array, 
                                            const uint32_t z_stagesPerStep, uint32_t z_rebasis)
{
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
    for(uint32_t i = 1; i < sim->n; i++)
    {  
      // Calculate the dt value and wrap around the orbital period.
      dt = fmod(dt, p_uVars->period[i]);

      // Calculate our step since last time we were called and update storage of tLast.
      double h = t - p_uVars->tLast[i];
      p_uVars->tLast[i] = t;

      SolveForUnivsersalAnomaly(dt, h, i, C);

      p_uVars->C.c0[i] = C[0];
      p_uVars->C.c1[i] = C[1];
      p_uVars->C.c2[i] = C[2];
      p_uVars->C.c3[i] = C[3];
    }

    for(uint32_t i = 1; i < sim->n; i++)
    {    
      double * Qout = Xosc_map[stage];
      double * Pout = &Qout[3*sim->n];      
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

      const double mi =  p_uVars->mass[i]; 
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
          RebasisOsculatingOrbits_Momenta(Qout, Pout, t, i);
        }
      }
    }
    
    if(z_rebasis != 0)
    {
      t_last_rebasis = h_array[stage];
    }
  } 
}


void ApplyCorrectorToOsculatingOrbitCalculation(double **Xosc_map, double t, uint32_t z_stagePerStep)
{ 
    double * Qout = Xosc_map[z_stagePerStep-1];
    double * Pout = &Qout[3*sim->n];      

    for(uint32_t i = 1; i < sim->n; i++)
    {
        for(uint32_t j = 0; j < 3; j++)
        {
          double csq = 0;
          double csp = 0;

          // Correct the osculating orbit inline with the various cs vars.
          add_cs(&Qout[i*3+j], &csq, -sim->radau->cs_dq[3*i+j]);
          add_cs(&Qout[i*3+j], &csq, -p_uVars->uv_csq[3*i+j]);
          add_cs(&Pout[i*3+j], &csp, -sim->radau->cs_dp[3*i+j]);
          add_cs(&Pout[i*3+j], &csp, -p_uVars->uv_csp[3*i+j]);   

          p_uVars->uv_csq[3*i+j] = 0;
          p_uVars->uv_csp[3*i+j] = 0;  
          sim->radau->cs_dq[3*i+j] = 0;
          sim->radau->cs_dp[3*i+j] = 0;        
          
          // Take anything left over from the oscualting orbit corrector and place into the delta.
          add_cs(&sim->radau->dQ[3*i+j], &p_uVars->uv_csq[3*i+j], -csq);
          add_cs(&sim->radau->dP[3*i+j], &p_uVars->uv_csp[3*i+j], -csp);

          // Obtain cs sum for v for next step.
          p_uVars->uv_csv[3*i+j] = p_uVars->uv_csp[3*i+j]/sim->mass[i];
        }
        RebasisOsculatingOrbits_Momenta(Qout, Pout, t, i);  
    }
  CalculateClassicalOrbitalElements(); // Remove this (All classical elements can be removed actually)
}

/*
This version is called whenever a rectification is performed.
*/
void RebasisOsculatingOrbits_Momenta(double * z_Q, double * z_P, double z_t, uint32_t i)
{

  p_uVars->Q0[3*i+0] = z_Q[3*i+0];
  p_uVars->Q0[3*i+1] = z_Q[3*i+1];
  p_uVars->Q0[3*i+2] = z_Q[3*i+2];
  
  p_uVars->V0[3*i+0] = z_P[3*i+0] / sim->mass[i];
  p_uVars->V0[3*i+1] = z_P[3*i+1] / sim->mass[i];
  p_uVars->V0[3*i+2] = z_P[3*i+2] / sim->mass[i];

  p_uVars->uv_csv[3*i] = p_uVars->uv_csp[3*i] / sim->mass[i];
  p_uVars->uv_csv[3*i+1] = p_uVars->uv_csp[3*i+1] / sim->mass[i];
  p_uVars->uv_csv[3*i+2] = p_uVars->uv_csp[3*i+2] / sim->mass[i];

  p_uVars->P0[3*i+0] = z_P[3*i+0];
  p_uVars->P0[3*i+1] = z_P[3*i+1];
  p_uVars->P0[3*i+2] = z_P[3*i+2];

  const double mu = p_uVars->mu;
  const double PIx4_MU = PI_SQ_X4 / mu;

  const double Q0x = z_Q[3*i];
  const double Q0y = z_Q[3*i+1];
  const double Q0z = z_Q[3*i+2];
  const double V0x = z_P[3*i]/sim->mass[i];
  const double V0y = z_P[3*i+1]/sim->mass[i];
  const double V0z = z_P[3*i+2]/sim->mass[i];

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

  CalculateClassicalOrbitalElementsSingle(i);
}


static double SolveForUnivsersalAnomaly(double dt, double h, uint32_t i, double * C)
{
  double X = p_uVars->X[i] + (h / p_uVars->Q0_norm[i]);
  X *= (1.0 - X*p_uVars->eta[i]*0.5/p_uVars->Q0_norm[i]);

  X = fmod(X, p_uVars->Xperiod[i]);

  double prevVals[MAX_NEWTON_ITERATIONS];

  for(uint32_t k = 0; k < MAX_NEWTON_ITERATIONS; k++)
  {    
    X = CalcUpdateValue(X, dt, i, C);

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

static double CalcUpdateValue(double X, double dt, uint32_t i, double * C)
{
  const double X2 = X*X;
  const double X3 = X2*X;
  double Z = p_uVars->beta[i]*X2;

  C_Stumpff(C, Z);

  const double num = (X*(p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2]) -
              p_uVars->eta[i]*X2*C[2] -
              p_uVars->zeta[i]*X3*C[3] + dt);
  const double denom = p_uVars->Q0_norm[i] + p_uVars->eta[i]*X*C[1] + p_uVars->zeta[i]*X2*C[2];

  return (num / denom);
}



static void C_Stumpff(double * cs, double z)
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

static void CalculateClassicalOrbitalElementsSingle(uint32_t i)
{
  if(sim->termination_check_enable)
  {
    // Semimajor axis
    p_uVars->a[i] = p_uVars->mu / p_uVars->beta[i];

    // Angular momentum
    CROSS(p_uVars->h, p_uVars->Q0, p_uVars->V0, i);
    p_uVars->h_norm[i] = NORM(p_uVars->h, i);

    // Eccentricity
    double eccVec[3] = {0,0,0};
    CROSS_LOCAL(eccVec, &p_uVars->V0[3*i], &p_uVars->h[3*i]);
    eccVec[0] /= p_uVars->mu;
    eccVec[1] /= p_uVars->mu;
    eccVec[2] /= p_uVars->mu;
    double Q_norm = NORM(p_uVars->Q0, i);
    eccVec[0] -= p_uVars->Q0[3*i+0] / Q_norm;
    eccVec[1] -= p_uVars->Q0[3*i+1] / Q_norm;
    eccVec[2] -= p_uVars->Q0[3*i+2] / Q_norm;
    p_uVars->e[i] = NORM(eccVec, 0);

    // Apogee and perigee.
    p_uVars->peri[i] = p_uVars->a[i]*(1-p_uVars->e[i]);
    p_uVars->apo[i] = p_uVars->a[i]*(1+p_uVars->e[i]);
  }
}

void CalculateClassicalOrbitalElements(void)
{
  if(sim->termination_check_enable)
  {
    for(uint32_t i = 1; i < sim->n; i++)
    {
      CalculateClassicalOrbitalElementsSingle(i);
    }
  }
}




void UniversalVars_Init(SIMULATION * z_sim)
{
  sim = z_sim;
  // Create the main control data structure for universal variables.
  p_uVars = (UNIVERSAL_VARS *)malloc(sizeof(UNIVERSAL_VARS));
  memset(p_uVars, 0, sizeof(UNIVERSAL_VARS));

  z_sim->uVars = p_uVars;

  p_uVars->stateVectorSize = 3 * z_sim->n * sizeof(double);
  p_uVars->controlVectorSize = z_sim->n * sizeof(double);

  // Allocate memory for all objects used in universal vars
  p_uVars->Q0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->V0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->Q1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->V1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->P0 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->P1 = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->t0 = (double *)malloc(z_sim->n*sizeof(double));
  p_uVars->tLast = (double *)malloc(z_sim->n*sizeof(double));
  p_uVars->Q0_norm = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->beta = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->eta =  (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->zeta = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->period = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->Xperiod = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->X = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->dt = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->mass = (double *)malloc(p_uVars->controlVectorSize);

  p_uVars->mu = (double)sim->G*(double)sim->mass[0];
  p_uVars->e = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->a = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->h = (double *)malloc(p_uVars->stateVectorSize);
  p_uVars->h_norm = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->peri = (double *)malloc(p_uVars->controlVectorSize);
  p_uVars->apo = (double *)malloc(p_uVars->controlVectorSize);
  
  // @todo dont need this copy - update all refs 
  for(uint32_t i = 0; i < z_sim->n; i++)
  {
    p_uVars->mass[i] = z_sim->mass[i];
  }

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
  memset(p_uVars->t0, 0, z_sim->n*sizeof(double));
  memset(p_uVars->tLast, 0, z_sim->n*sizeof(double));
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
  p_uVars->uv_csq = (double *)malloc(3 * z_sim->n * sizeof(double));
  p_uVars->uv_csv = (double *)malloc(3 * z_sim->n * sizeof(double));
  p_uVars->uv_csp = (double *)malloc(3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csq, 0, 3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csv, 0, 3 * z_sim->n * sizeof(double));
  memset(p_uVars->uv_csp, 0, 3 * z_sim->n * sizeof(double));  
}


void UniversalVars_Free(void)
{
  free(p_uVars->Q0);
  free(p_uVars->V0);
  free(p_uVars->Q1);
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
  free(p_uVars);
  sim->uVars = NULL;
}

static inline void add_cs(double * out, double * cs, double inp)
{
    const double y = inp - cs[0];
    const double t = out[0] + y;
    cs[0] = (t - out[0]) - y;
    out[0] = t;      
}
