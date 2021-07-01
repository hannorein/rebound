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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "radau.h"
#include "radau_step.h"
#include <math.h>
#include "dhem.h"
#include "UniversalVars.h"

static SIMULATION * __restrict__ sim;
static RADAU * __restrict__ radau;

// Note this is modified from previous implementations and now has 1.0 instead of 0.0.
double hArr[9] = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1.0};
static const double rr_inv[28] = {17.7738089140780008407526624, 5.5481367185372165056928203, 8.0659386483818866885371223, 2.8358760786444386782520107, 3.3742499769626352599420361, 5.8010015592640614823286804, 1.8276402675175978297946079, 2.0371118353585847827949161, 2.7254422118082262837742729, 5.1406241058109342286363203, 1.3620078160624694969370006, 1.4750402175604115479218482, 1.8051535801402512604391149, 2.6206449263870350811541814, 5.3459768998711075141214895, 1.1295338753367899027322862, 1.2061876660584456166252037, 1.4182782637347391537713785, 1.8772424961868100972169920, 2.9571160172904557478071039, 6.6176620137024244874471300, 1.0229963298234867458386119, 1.0854721939386423840467243, 1.2542646222818777659905423, 1.6002665494908162609916716, 2.3235983002196942228325344, 4.1099757783445590862385765, 10.8460261902368446847064289};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};
static const double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588};


#define stages 7
#define OSCULATING_ORBIT_SLOTS 9  // stages + 2 for t=0 and t=1.
#define MAX_ITERATIONS 12
#define MIN_ITERATIONS 2
#define NORM(x, i) sqrt(x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2])

#define LOCAL_ERROR_ESTIMATE 0
#define GLOBAL_ERROR_ESTIMATE 1
// #define ERROR_ESTIMATE LOCAL_ERROR_ESTIMATE
#define ERROR_ESTIMATE GLOBAL_ERROR_ESTIMATE

#define CUT_OFF 1E-13

#define DONT_ADD_CS_VAR 0
#define ADD_CS_VAR 1

// Integration step variables
controlVars G;
controlVars B;
controlVars Blast;

controlVars B0;
controlVars Blast0;
controlVars Blast_1st;

controlVars * bActive;

controlVars G_1st;
controlVars B_1st;

controlVars G_full;
controlVars B_full;
controlVars Blast_full;


static void RadauStep15_Step(uint32_t * z_fCalls, double t, double h, uint32_t step);
static inline void add_cs(double* out, double* cs, double inp);
static double ReturnIAS15StepError(double h, double t);
static void ControlVars_Copy(controlVars * out, controlVars * in);
static void CalculateGFromBInternal(controlVars * z_G, controlVars * z_B, uint32_t z_start, uint32_t z_end);


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

controlVars_const control_vars_cast(controlVars in)
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

static void RadauStep15_Step(uint32_t * z_iterations, double t, double h, uint32_t step)
{
  double errMax = 0;
  controlVars * z_csB = &radau->cs_B;
  const uint32_t n3 = sim->n3;

  sim->f_rhs(radau->dQ, radau->dP, radau->dState0, &radau->dState0[3*sim->n],
             radau->ddState0, &radau->ddState0[3*sim->n], 0, radau->cs_dState0, radau->cs_ddState0);

  radau->fCalls++;

  // Wipe the cs variables for B summation.
  ControlVars_Clear(z_csB);
  ControlVars_Clear(&radau->cs_B1st);

  uint32_t n = 1;;
  for(n = 1; n < MAX_ITERATIONS; n++)
  {
    z_iterations[0] = n + 1;

    for(uint32_t i = 1; i <= stages; i++)
    {
      radau->t = t+h*hArr[i];
      CalculatePredictors(h, hArr[i], radau->dX, radau->dState0, radau->ddState0, &B,
                          radau->predictors, radau->cs_dX, 3, sim->stateVectorLength/2, DONT_ADD_CS_VAR);

      CalculatePredictors_1stOrder(h, hArr[i], radau->dX, radau->dState0, &B_1st, radau->predictors, 
                                  radau->cs_dX, (int)sim->stateVectorLength/2, sim->stateVectorLength, DONT_ADD_CS_VAR);


      sim->f_rhs(radau->predictors, &radau->predictors[3*sim->n], radau->dState, &radau->dState[3*sim->n],
           radau->ddState, &radau->ddState[3*sim->n], i,
           radau->cs_dState, radau->cs_ddState);

      radau->fCalls++;

      switch(i)
      {
        case 1: 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p0[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];                  
                  G.p0[i] = Gi*rr_inv[0];

                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), G.p0[i]-tmp); 
                }                

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p0[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st.p0[j] = Gi2*rr_inv[0];
                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), G_1st.p0[j]-tmp2);
                }                       
                break;
        case 2: 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p1[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G.p1[i] = (Gi*rr_inv[1] - G.p0[i]) *rr_inv[2];
                  tmp = G.p1[i] - tmp;
                  
                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[0]);
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp);
                }    

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p1[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st.p1[j] = (Gi2*rr_inv[1] - G_1st.p0[j]) *rr_inv[2];
                  tmp2 = G_1st.p1[j] - tmp2;

                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[0]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2);                    
                }                
                break;
        case 3:                 
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p2[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  G.p2[i] = ((Gi*rr_inv[3] - G.p0[i]) *rr_inv[4] - G.p1[i])*rr_inv[5];

                  tmp = G.p2[i] - tmp;

                  // We only need to calculate the change in B.
                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[1]);
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp*c[2]);
                  add_cs(&(B.p2[i]), &(radau->cs_B.p2[i]), tmp);            
                }       

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p2[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st.p2[j] = ((Gi2*rr_inv[3] - G_1st.p0[j]) *rr_inv[4] - G_1st.p1[j])*rr_inv[5];

                  tmp2 = G_1st.p2[j] - tmp2;

                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[1]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[2]);
                  add_cs(&(B_1st.p2[j]), &(radau->cs_B1st.p2[j]), tmp2);      
                }                         
                break;
        case 4:  
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p3[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];

                  G.p3[i] = (((Gi*rr_inv[6] - G.p0[i]) *rr_inv[7] - G.p1[i])*rr_inv[8] - G.p2[i]) *rr_inv[9];

                  tmp = G.p3[i] - tmp;

                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[3]);
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp*c[4]);
                  add_cs(&(B.p2[i]), &(radau->cs_B.p2[i]), tmp*c[5]);
                  add_cs(&(B.p3[i]), &(radau->cs_B.p3[i]), tmp);              
                }     

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p3[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];

                  G_1st.p3[j] = (((Gi2*rr_inv[6] - G_1st.p0[j]) *rr_inv[7] - G_1st.p1[j])*rr_inv[8] - G_1st.p2[j]) *rr_inv[9];

                  tmp2 = G_1st.p3[j] - tmp2;

                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[3]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[4]);
                  add_cs(&(B_1st.p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[5]);
                  add_cs(&(B_1st.p3[j]), &(radau->cs_B1st.p3[j]), tmp2);    
                }                

                break;
        case 5:  
                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p4[i];
                  double Gi = radau->ddState[i]-radau->ddState0[i];                  
                  G.p4[i] = ((((Gi*rr_inv[10] - G.p0[i]) *rr_inv[11] - G.p1[i])*rr_inv[12] - G.p2[i])*rr_inv[13] - G.p3[i]) *rr_inv[14];
                  tmp = G.p4[i] - tmp;
                  

                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[6]);
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp*c[7]);
                  add_cs(&(B.p2[i]), &(radau->cs_B.p2[i]), tmp*c[8]);
                  add_cs(&(B.p3[i]), &(radau->cs_B.p3[i]), tmp*c[9]);
                  add_cs(&(B.p4[i]), &(radau->cs_B.p4[i]), tmp);
                }

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p4[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];
                  G_1st.p4[j] = ((((Gi2*rr_inv[10] - G_1st.p0[j]) *rr_inv[11] - G_1st.p1[j])*rr_inv[12] - G_1st.p2[j])*rr_inv[13] - G_1st.p3[j]) *rr_inv[14];
                  tmp2 = G_1st.p4[j] - tmp2;
                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[6]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[7]);
                  add_cs(&(B_1st.p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[8]);
                  add_cs(&(B_1st.p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[9]);
                  add_cs(&(B_1st.p4[j]), &(radau->cs_B1st.p4[j]), tmp2);          
                }                

                break;
        case 6:  
                #pragma GCC ivdep 
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p5[i];

                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G.p5[i] = (((((Gi*rr_inv[15] - G.p0[i]) *rr_inv[16] - G.p1[i])*rr_inv[17] -
                                    G.p2[i]) *rr_inv[18] - G.p3[i]) *rr_inv[19] - G.p4[i]) *rr_inv[20];
                  
                  tmp = G.p5[i] - tmp;
                  

                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[10]);                  
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp*c[11]);
                  add_cs(&(B.p2[i]), &(radau->cs_B.p2[i]), tmp*c[12]);
                  add_cs(&(B.p3[i]), &(radau->cs_B.p3[i]), tmp*c[13]);
                  add_cs(&(B.p4[i]), &(radau->cs_B.p4[i]), tmp*c[14]);
                  add_cs(&(B.p5[i]), &(radau->cs_B.p5[i]), tmp);
                }

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p5[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j]; 
                  G_1st.p5[j] = (((((Gi2*rr_inv[15] - G_1st.p0[j]) *rr_inv[16] - G_1st.p1[j])*rr_inv[17] -
                                    G_1st.p2[j]) *rr_inv[18] - G_1st.p3[j]) *rr_inv[19] - G_1st.p4[j]) *rr_inv[20]; 
                  tmp2 = G_1st.p5[j] - tmp2;

                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[10]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[11]);
                  add_cs(&(B_1st.p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[12]);
                  add_cs(&(B_1st.p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[13]);
                  add_cs(&(B_1st.p4[j]), &(radau->cs_B1st.p4[j]), tmp2*c[14]);
                  add_cs(&(B_1st.p5[j]), &(radau->cs_B1st.p5[j]), tmp2);

                }
                break;
        case 7:  
                radau->acc_ptr = radau->ddState;
                bActive = &B;

                #pragma GCC ivdep
                for(uint32_t j = n3+3; j < 2*n3; j++)
                {
                  double tmp2 = G_1st.p6[j];
                  double Gi2 = radau->dState[j]-radau->dState0[j];

                  G_1st.p6[j] = ((((((Gi2*rr_inv[21] - G_1st.p0[j]) *rr_inv[22] - G_1st.p1[j])*rr_inv[23] - G_1st.p2[j]) *rr_inv[24] - 
                                  G_1st.p3[j]) *rr_inv[25] - G_1st.p4[j]) *rr_inv[26] - G_1st.p5[j]) *rr_inv[27];
                  tmp2 = G_1st.p6[j] - tmp2;                                  

                  add_cs(&(B_1st.p0[j]), &(radau->cs_B1st.p0[j]), tmp2*c[15]);
                  add_cs(&(B_1st.p1[j]), &(radau->cs_B1st.p1[j]), tmp2*c[16]);
                  add_cs(&(B_1st.p2[j]), &(radau->cs_B1st.p2[j]), tmp2*c[17]);
                  add_cs(&(B_1st.p3[j]), &(radau->cs_B1st.p3[j]), tmp2*c[18]);
                  add_cs(&(B_1st.p4[j]), &(radau->cs_B1st.p4[j]), tmp2*c[19]);
                  add_cs(&(B_1st.p5[j]), &(radau->cs_B1st.p5[j]), tmp2*c[20]);
                  add_cs(&(B_1st.p6[j]), &(radau->cs_B1st.p6[j]), tmp2);
                }

                #pragma GCC ivdep
                for(uint32_t i = 3; i < n3; i++)
                {
                  double tmp = G.p6[i]; 
                  double Gi = radau->ddState[i]-radau->ddState0[i];
                  
                  G.p6[i] =     ((((((Gi*rr_inv[21] - G.p0[i]) *rr_inv[22] - G.p1[i])*rr_inv[23] - G.p2[i]) *rr_inv[24] - 
                                  G.p3[i]) *rr_inv[25] - G.p4[i]) *rr_inv[26] - G.p5[i]) *rr_inv[27];
                                  
                  tmp = G.p6[i] - tmp;
                  
                  add_cs(&(B.p0[i]), &(radau->cs_B.p0[i]), tmp*c[15]);
                  add_cs(&(B.p1[i]), &(radau->cs_B.p1[i]), tmp*c[16]);
                  add_cs(&(B.p2[i]), &(radau->cs_B.p2[i]), tmp*c[17]);
                  add_cs(&(B.p3[i]), &(radau->cs_B.p3[i]), tmp*c[18]);
                  add_cs(&(B.p4[i]), &(radau->cs_B.p4[i]), tmp*c[19]);
                  add_cs(&(B.p5[i]), &(radau->cs_B.p5[i]), tmp*c[20]);
                  add_cs(&(B.p6[i]), &(radau->cs_B.p6[i]), tmp);
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
    double q_ddot[3*sim->n];

    for(uint32_t j = 1; j < sim->n; j++)
    {
      q_ddot[3*j+0] = sim->rhs->Xosc_dotArr[8][3*sim->n+3*j+0] / sim->mass[j];
      q_ddot[3*j+1] = sim->rhs->Xosc_dotArr[8][3*sim->n+3*j+1] / sim->mass[j];
      q_ddot[3*j+2] = sim->rhs->Xosc_dotArr[8][3*sim->n+3*j+2] / sim->mass[j];
    }    

    for(uint32_t j = 3; j < 3*sim->n; j++)
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
  

    for(uint32_t k = 1; k < sim->n; k++)
    {
      ratio = NORM(radau->dX, k) / NORM(sim->rhs->Qosc, k);

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
  CalculateNewState(h, radau->dState0, radau->ddState0, &B, radau->dX, 
                    radau->cs_dX, 0, (int)(n3));

  CalculateNewState_1stOrder(h, radau->dState0, &B_1st, radau->dX, 
                                radau->cs_dX, (int)(n3), sim->stateVectorLength);                                

  // Now that we have cs_dX we can perform our correction.
  ApplyCorrectorToOsculatingOrbitCalculation(sim->rhs->XoscArr, t+h, 9); 
}

static double ReturnIAS15StepError(double h, double t)
{
  double b6Max = 0;
  double accMax = 0;
  double errMax = 0;

  for(uint32_t i = 0; i < sim->n; i++) 
  {
      for(uint32_t j = 0; j < 3; j++)
      {
        double b6 = fabs(bActive->p6[3*i+j]);
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

  radau->b6Max = b6Max;
  radau->accMax = accMax;

  return errMax;
}

static void CalculateGFromBInternal(controlVars * z_G, controlVars * z_B, uint32_t z_start, uint32_t z_end)
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

void CalculateGfromB(void)
{
  CalculateGFromBInternal(&G, &B, 0, (int)(sim->stateVectorLength/2));
  CalculateGFromBInternal(&G_1st, &B_1st, (int)(sim->stateVectorLength/2), sim->stateVectorLength);
}

void AnalyticalContinuation(controlVars * z_B, controlVars * z_Blast, const double h,
                             const double h_new, const uint32_t * const rectificationArray, const uint32_t step)
{
  const double ratio = h_new / h;
  const double q1 = ratio;
  const double q2 = q1 * q1;
  const double q3 = q1 * q2;
  const double q4 = q2 * q2;
  const double q5 = q2 * q3;
  const double q6 = q3 * q3;
  const double q7 = q3 * q4;

  for(uint32_t i = 0; i < sim->stateVectorLength; i++)
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

void CalculatePredictors_1stOrder(double h, double hSample, double * z_state0, double * z_dState0, 
                                    controlVars * z_B, double * z_predictors, double * z_csState, 
                                    uint32_t z_start, uint32_t z_end, uint32_t add_cs_var)
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

            
void CalculatePredictors(double h, double hSample, double const * __restrict__ z_state0, double const * __restrict__ z_dState, 
                         double const * __restrict__ z_ddState, controlVars const * z_B, double * __restrict__ z_predictors, 
                         double const * __restrict__ z_csState, uint32_t const z_start, uint32_t const z_end, uint32_t const add_cs_var)
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






void CalculateNewState_1stOrder(double h, double * z_dState, controlVars * z_B, double * z_state0, 
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

void CalculateNewState(double h, double * z_dState, double * z_ddState,
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

void RejectStep(void)
{
  memcpy(radau->dX, radau->dX0, sim->stateVectorSize);
  ControlVars_Copy(&B, &B0);
  ControlVars_Copy(&Blast, &Blast0);
}

void RadauStep15_Init(SIMULATION * z_sim)
{
    sim = z_sim;
    radau = sim->radau;
    radau->step = RadauStep15_Step;
    radau->AnalyticalContinuation = AnalyticalContinuation;
    radau->CalculateGfromB = CalculateGfromB;
    radau->ReturnStepError = ReturnIAS15StepError;
    radau->RejectStep = RejectStep;

    ControlVars_Init(&G, sim->stateVectorSize);
    ControlVars_Init(&B, sim->stateVectorSize);
    ControlVars_Init(&radau->cs_B, sim->stateVectorSize);
    ControlVars_Init(&Blast, sim->stateVectorSize);
    ControlVars_Init(&Blast_1st, sim->stateVectorSize);

    ControlVars_Init(&B0, sim->stateVectorSize);
    ControlVars_Init(&Blast0, sim->stateVectorSize);
    ControlVars_Init(&G_full, sim->stateVectorSize);
    ControlVars_Init(&B_full, sim->stateVectorSize);
    ControlVars_Init(&Blast_full, sim->stateVectorSize);

    radau->B = &B;
    radau->Blast = &Blast;

    radau->B_1st = &B_1st;
    radau->Blast_1st = &Blast_1st;

    radau->B_full = &B_full;
    radau->Blast_full = &Blast_full;

    ControlVars_Init(&G_1st, sim->stateVectorSize);
    ControlVars_Init(&B_1st, sim->stateVectorSize);
    ControlVars_Init(&radau->cs_B1st, sim->stateVectorSize);
    ControlVars_Init(&radau->cs_B_full, sim->stateVectorSize);

    radau->dState0 = (double *)malloc(sim->stateVectorSize);
    radau->ddState0 = (double *)malloc(sim->stateVectorSize);
    radau->dState = (double *)malloc(sim->stateVectorSize);
    radau->ddState = (double *)malloc(sim->stateVectorSize);

    memset(radau->dState0, 0, sim->stateVectorSize);
    memset(radau->ddState0, 0, sim->stateVectorSize);
    memset(radau->dState, 0, sim->stateVectorSize);
    memset(radau->ddState, 0, sim->stateVectorSize);

    // Compensated summation arrays
    radau->cs_state = (double *)malloc(sim->stateVectorSize);
    radau->csx = (double *)malloc(sim->stateVectorSize);
    radau->csv = (double *)malloc(sim->stateVectorSize);
    radau->cs_state_1st = (double *)malloc(sim->stateVectorSize);
    radau->cs_dState0 = (double *)malloc(sim->stateVectorSize);
    radau->cs_ddState0 = (double *)malloc(sim->stateVectorSize);
    radau->cs_dState = (double *)malloc(sim->stateVectorSize);
    radau->cs_ddState = (double *)malloc(sim->stateVectorSize);

    memset(radau->cs_state, 0, sim->stateVectorSize);
    memset(radau->csx, 0, sim->stateVectorSize);
    memset(radau->csv, 0, sim->stateVectorSize);
    memset(radau->cs_state_1st, 0, sim->stateVectorSize);
    memset(radau->cs_dState0, 0, sim->stateVectorSize);
    memset(radau->cs_ddState0, 0, sim->stateVectorSize);
    memset(radau->cs_dState, 0, sim->stateVectorSize);
    memset(radau->cs_ddState, 0, sim->stateVectorSize);


    radau->cs_dq_dot = (double *)malloc(sim->stateVectorSize);
    radau->cs_dX = (double *)malloc(sim->stateVectorSize);
    radau->cs_dq = radau->cs_dX;
    radau->cs_dp = &radau->cs_dX[(int)sim->stateVectorLength/2];

    radau->csx_recti = (double *)malloc(sim->stateVectorSize);
    radau->csv_recti = (double *)malloc(sim->stateVectorSize);
    
    memset(radau->cs_dq_dot, 0, sim->stateVectorSize);
    memset(radau->cs_dX, 0, sim->stateVectorSize);
    memset(radau->csx_recti, 0, sim->stateVectorSize);
    memset(radau->csv_recti, 0, sim->stateVectorSize);
}

void RadauStep15_Free(void)
{
  ControlVars_Free(&G);
  ControlVars_Free(&B);
  ControlVars_Free(&radau->cs_B);
  ControlVars_Free(&Blast);
  ControlVars_Free(&B0);
  ControlVars_Free(&Blast0);
  free(radau->dState0);
  free(radau->ddState0);
  free(radau->dState);
  free(radau->ddState);
  free(radau->cs_dState0);
  free(radau->cs_ddState0);
  free(radau->cs_dState);
  free(radau->cs_ddState);
  free(radau->cs_state);
  free(radau->cs_state_1st);
  // free(radau->dState0_hp);
}

static void ControlVars_Copy(controlVars * out, controlVars * in)
{
  if(out->size == in->size)
  {
    memcpy(out->p0, in->p0, out->size);
    memcpy(out->p1, in->p1, out->size);
    memcpy(out->p2, in->p2, out->size);
    memcpy(out->p3, in->p3, out->size);
    memcpy(out->p4, in->p4, out->size);
    memcpy(out->p5, in->p5, out->size);
    memcpy(out->p6, in->p6, out->size);
  }
}

void ControlVars_Init(controlVars * var, uint32_t size)
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

void ControlVars_Free(controlVars * var)
{
  free(var->p0);
  free(var->p1);
  free(var->p2);
  free(var->p3);
  free(var->p4);
  free(var->p5);
  free(var->p6);
}

void ControlVars_Clear(controlVars * var)
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
  var->size = 0;
}


static inline void add_cs(double* out, double* cs, double inp)
{
      const double y = inp - cs[0];
      const double t = out[0] + y;
      cs[0] = (t - out[0]) - y;
      out[0] = t;      
}
