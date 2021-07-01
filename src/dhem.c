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

#include "radau.h"
#include "dhem.h"
#include "UniversalVars.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static DHEM * dhem;
static SIMULATION * sim;

#define NORM(x, i) sqrt(x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2])
#define NORMl(x, i) sqrtq(x[3*i]*x[3*i]+x[3*i+1]*x[3*i+1]+x[3*i+2]*x[3*i+2])
#define DOT(x, y, i, j) (x[3*i+0]*y[3*j+0]+x[3*i+1]*y[3*j+1]+x[3*i+2]*y[3*j+2])

static void dhem_PerformSummation(double * Q, double * P,
                               double * dQ, double * dP, uint32_t stageNumber);
static void dhem_InitialiseOsculatingOrbits(double * Q, double * P, double t);
static inline void add_cs(double* out, double* cs, double inp);

void dhem_rhs_wrapped(double * dQ, double * dP, double * dQ_dot,
              double * dP_dot, double * dQ_ddot, double * dP_ddot, uint32_t stageNumber,
              double * cs1, double * cs2)
{
  // Set up the pointer to the previously calculated osculating orbit values.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*sim->n];
  dhem->Xosc_dot = dhem->Xosc_dotArr[stageNumber];
  dhem->Qosc_dot = dhem->Xosc_dot;
  dhem->Posc_dot = &dhem->Qosc_dot[3*sim->n];

  // cs vars
  dhem->Xosc_cs = dhem->XoscArr_cs[stageNumber];
  dhem->Qosc_cs = dhem->Xosc_cs;
  dhem->Posc_cs = &dhem->Qosc_cs[3*sim->n];

  dhem_rhs(dQ, dP, dQ_dot, dP_dot, dQ_ddot, dP_ddot);

}


/*
 * dQ, dP are the input deltas and dQ_dot, dP_dot are the first derivatives terms.
 * dQ_ddot and dP_ddot are the second derivative terms.
 */
void dhem_rhs(double const * __restrict__ const dQ, double const * __restrict__ const dP, double * __restrict__ const dQ_dot,
              double * __restrict__ const dP_dot, double * __restrict__ const dQ_ddot, double * __restrict__ const dP_ddot)
{
  // Not necessary but makes code more reable.
  const double * const __restrict__ m = dhem->m;
  const uint32_t n = dhem->n;
  const uint32_t n3 = 3*dhem->n;
  const double G = sim->G;
  const double * const __restrict__ Qosc = dhem->Qosc;
  const double * const __restrict__ Posc = dhem->Posc;
  const double * const __restrict__ Posc_dot = dhem->Posc_dot;
  double * __restrict__ Q = dhem->Q;
  double * __restrict__ P = dhem->P;

  memset(dQ_ddot, 0, sim->stateVectorSize/2);

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

double dhem_CalculateHamiltonian(double * Q, double * P)
{
  double * m = sim->mass;
  double hamiltonian = 0;
  double Psun[3] = {0,0,0};
  double sepVec[3] = {0,0,0};

  for(uint32_t i = 1; i < sim->n; i++)
  {
    double Pnorm = (NORM(P, i));
    hamiltonian += Pnorm*Pnorm / (2*m[i]);
    hamiltonian -= sim->G*m[0]*m[i] / NORM(Q, i);

    Psun[0] += P[3*i];
    Psun[1] += P[3*i+1];
    Psun[2] += P[3*i+2];
  }
  double PsunNorm = NORM(Psun, 0);
  hamiltonian += (1.0/(2.0*m[0])) * PsunNorm*PsunNorm;

  for(uint32_t i = 1; i < sim->n; i++)
  {
    for(uint32_t j = i+1; j < sim->n; j++)
    {
      sepVec[0] = Q[3*i]-Q[3*j];
      sepVec[1] = Q[3*i+1]-Q[3*j+1];
      sepVec[2] = Q[3*i+2]-Q[3*j+2];

      hamiltonian -= (sim->G*m[i]*m[j]) / NORM(sepVec, 0);
    }
  }

  return hamiltonian;
}

static void dhem_PerformSummation(double * Q, double * P,
                               double * dQ, double * dP, uint32_t stageNumber)
{
  // Need to setup a pointer to the osculating orbits in memory.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*sim->n];

  // Always zero in our frame of reference.
  for(uint32_t i = 0; i < 3; i++)
  {
    Q[i] = 0;
    P[i] = 0;
  }

  for(uint32_t i = 1; i < sim->n; i++)
  {
    Q[3*i+0] = dhem->Qosc[3*i+0] + (dQ[3*i+0] + (dhem->Qosc_cs[3*i+0] + sim->radau->cs_dq[3*i+0]));
    Q[3*i+1] = dhem->Qosc[3*i+1] + (dQ[3*i+1] + (dhem->Qosc_cs[3*i+1] + sim->radau->cs_dq[3*i+1]));
    Q[3*i+2] = dhem->Qosc[3*i+2] + (dQ[3*i+2] + (dhem->Qosc_cs[3*i+2] + sim->radau->cs_dq[3*i+2]));

    P[3*i+0] = dhem->Posc[3*i+0] + (dP[3*i+0] + (dhem->Posc_cs[3*i+0] + sim->radau->cs_dp[3*i+0]));
    P[3*i+1] = dhem->Posc[3*i+1] + (dP[3*i+1] + (dhem->Posc_cs[3*i+1] + sim->radau->cs_dp[3*i+1]));
    P[3*i+2] = dhem->Posc[3*i+2] + (dP[3*i+2] + (dhem->Posc_cs[3*i+2] + sim->radau->cs_dp[3*i+2]));
  }

}

static void dhem_InitialiseOsculatingOrbits(double * Q, double * P, double t)
{
  for(uint32_t i = 1; i < sim->n; i++)
  {
    RebasisOsculatingOrbits_Momenta(Q, P, t, i);
  }
}

uint32_t dhem_RectifyOrbits(double t, double * Q, double * P,
                            double * dQ, double * dP, uint32_t * rectifiedArray, uint32_t stageNumber)
{
  uint32_t rectifyFlag = 0;
  uint32_t rectifiedCount = 0;
  double dQ_norm = 0;
  double dP_norm = 0;

  // Need to setup a pointer to the osculating orbits in memory.
  dhem->Xosc = dhem->XoscArr[stageNumber];
  dhem->Qosc = dhem->Xosc;
  dhem->Posc = &dhem->Qosc[3*sim->n];
  dhem->Xosc_dot = dhem->Xosc_dotArr[stageNumber];
  dhem->Qosc_dot = dhem->Xosc_dot;
  dhem->Posc_dot = &dhem->Qosc_dot[3*sim->n];

  // CS variables
  dhem->Xosc_cs = dhem->XoscArr_cs[stageNumber];
  dhem->Qosc_cs = dhem->Xosc_cs;
  dhem->Posc_cs = &dhem->Qosc_cs[3*sim->n];

  for(uint32_t i = 1; i < sim->n; i++)
  {
    rectifiedArray[3*i] = 0;
    rectifiedArray[3*i+1] = 0;
    rectifiedArray[3*i+2] = 0;

    rectifiedArray[sim->n*3+3*i] = 0;
    rectifiedArray[sim->n*3+3*i+1] = 0;
    rectifiedArray[sim->n*3+3*i+2] = 0;

    dQ_norm = NORM(dQ, i) / NORM(dhem->Qosc, i);
    dP_norm = NORM(dP, i) / NORM(dhem->Posc, i);

    if(t > dhem->rectifyTimeArray[i] ||
      dQ_norm > sim->dQcutoff ||
      dP_norm > sim->dPcutoff)
    {
      rectifyFlag = 1;
      break;
    }
  }

  for(uint32_t i = 1; i < sim->n; i++)
  {
    if(rectifyFlag != 0)
    {
      rectifiedCount++;

      for(uint32_t j = 0; j < 3; j++)
      {
        double temp_cs = 0;

        Q[3*i+j] = dhem->Qosc[3*i+j];
        add_cs(&Q[3*i+j], &temp_cs, dQ[3*i+j]);
        add_cs(&Q[3*i+j], &temp_cs, sim->uVars->uv_csq[3*i+j]);            
        add_cs(&Q[3*i+j], &temp_cs, sim->radau->cs_dq[3*i+j]);

        dQ[3*i+j] = -temp_cs;
        sim->radau->cs_dq[3*i+j] = 0;
        sim->uVars->uv_csq[3*i+j] = 0;
      }

      for(uint32_t j = 0; j < 3; j++)
      {
        double temp_cs = 0;

        P[3*i+j] = dhem->Posc[3*i+j];
        add_cs(&P[3*i+j], &temp_cs, dP[3*i+j]);
        add_cs(&P[3*i+j], &temp_cs, sim->uVars->uv_csp[3*i+j]);            
        add_cs(&P[3*i+j], &temp_cs, sim->radau->cs_dp[3*i+j]);

        dP[3*i+j] = -temp_cs;
        sim->radau->cs_dp[3*i+j] = 0;
        sim->uVars->uv_csp[3*i+j] = 0;
      }          

      RebasisOsculatingOrbits_Momenta(Q, P, t, i);
      //dhem->rectifyTimeArray[i] += dhem->rectificationPeriod[i];
      // This option will add some randomness to the rectification process.
      dhem->rectifyTimeArray[i] = t + dhem->rectificationPeriod[i];

      rectifiedArray[3*i] = 1;
      rectifiedArray[3*i+1] = 1;
      rectifiedArray[3*i+2] = 1;
      rectifiedArray[sim->n*3+3*i] = 1;
      rectifiedArray[sim->n*3+3*i+1] = 1;
      rectifiedArray[sim->n*3+3*i+2] = 1;
    }    
  }


  return rectifiedCount;
}

void dhem_CalculateOsculatingOrbitDerivatives_Momenta(double const * const __restrict__ Qosc, double const * const __restrict__ Posc, 
                                                      double * const __restrict__ Qosc_dot, double * const __restrict__ Posc_dot)
{
  const double GM0 = -sim->G*dhem->m[0];

  for(uint32_t i = 1; i < sim->n; i++)
  {
    const double m = sim->mass[i];
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


void dhem_ConvertToCOM(double * Q, double * V, double * Qout, double * Vout)
{
  double * m = sim->mass;
  double mTotal = 0;
  double COM[3] = {0,0,0};
  double COP[3] = {0,0,0};

  switch(sim->initial_condition_format)
  {
    case cartesian:
    {
        
      for(uint32_t i = 0; i < sim->n; i++)
      {
        mTotal += m[i];

        COM[0] += Q[3*i]*m[i];
        COM[1] += Q[3*i+1]*m[i];
        COM[2] += Q[3*i+2]*m[i];
      }
      COM[0] /= mTotal;
      COM[1] /= mTotal;
      COM[2] /= mTotal;

      for(uint32_t i = 0; i < sim->n; i++)
      {
        Qout[3*i] = Q[3*i] - COM[0];
        Qout[3*i+1] = Q[3*i+1] - COM[1];
        Qout[3*i+2] = Q[3*i+2] -COM[2];
      }

      // Centre of momentum
      for(uint32_t i = 0; i < sim->n; i++)
      {
        COP[0] += V[3*i]*m[i];
        COP[1] += V[3*i+1]*m[i];
        COP[2] += V[3*i+2]*m[i];
      }
      COP[0] /= mTotal;
      COP[1] /= mTotal;
      COP[2] /= mTotal;

      for(uint32_t i = 0; i < sim->n; i++)
      {
        Vout[3*i] = V[3*i] - COP[0]/m[i];
        Vout[3*i+1] = V[3*i+1] - COP[1]/m[i];
        Vout[3*i+2] = V[3*i+2] - COP[2]/m[i];
      }
      break;
    }
    case democratic_heliocentric:
      break;
  }

}

void dhem_ConvertToDHCoords(double * Q, double * V, double * Qout, double * Pout)
{
  double * m = sim->mass;

  switch(sim->initial_condition_format)
  {
    case cartesian:
    {
      printf("\ncartesian input.");
      // Calculate the Q terms other than Q0 as distance from the central body.
      for(uint32_t i = 1; i < sim->n; i++)
      {
        Qout[3*i] = Q[3*i] - Q[0];
        Qout[3*i+1] = Q[3*i+1] - Q[1];
        Qout[3*i+2] = Q[3*i+2] - Q[2];
      }

      double mTot = dhem->mTotal;

      // Work in the COM frame.
      Qout[0] = 0;
      Qout[1] = 0;
      Qout[2] = 0;

      // Calculate the momenta terms.
      double COP[3] = {0,0,0};

      for(uint32_t i = 0; i < sim->n; i++)
      {
        COP[0] += V[3*i]*m[i];
        COP[1] += V[3*i+1]*m[i];
        COP[2] += V[3*i+2]*m[i];
      }

      for(uint32_t i = 0; i < sim->n; i++)
      {
          V[3*i] -= COP[0] / mTot;
          V[3*i+1] -= COP[1] / mTot;
          V[3*i+2] -= COP[2] / mTot;
      }

      for(uint32_t i = 0; i < sim->n; i++)
      {
        Pout[3*i] = V[3*i]*m[i];
        Pout[3*i+1] = V[3*i+1]*m[i];
        Pout[3*i+2] = V[3*i+2]*m[i];
      }

      for(uint32_t i = 1; i < sim->n; i++)
      {
          Pout[0] += Pout[3*i];
          Pout[1] += Pout[3*i+1];
          Pout[2] += Pout[3*i+2];
      }

      for(uint32_t i = 1; i < sim->n; i++)
      {
        Pout[3*i] = Pout[3*i] - (m[i]/m[0])*Pout[0];
        Pout[3*i+1] = Pout[3*i+1] - (m[i]/m[0])*Pout[1];
        Pout[3*i+2] = Pout[3*i+2] - (m[i]/m[0])*Pout[2];
      }
      break;
    }
    case democratic_heliocentric:
    {
      printf("\ndemocratic heliocentric input.");
      for(uint32_t i = 0; i < 3*sim->n; i++)      
      {
        Qout[i] = Q[i];
        Pout[i] = V[i];
      }   
      break;
    }
  }
      // Work in the COM frame.
      Pout[0] = 0;
      Pout[1] = 0;
      Pout[2] = 0;
}

void dhem_CalcOscOrbitsForAllStages(double t0, double h, double * hArr, uint32_t z_stagesPerStep, uint32_t z_rebasis)
{
  if(z_rebasis != 0)
  {
    CalculateOsculatingOrbitsForSingleStep(dhem->XoscArr, t0, h, hArr, z_stagesPerStep, z_rebasis);  
  }
  else
  {
    CalculateOsculatingOrbitsForSingleStep(dhem->XoscPredArr, t0, h, hArr, z_stagesPerStep, z_rebasis);  
  }

  for(uint32_t i = 0; i < z_stagesPerStep; i++)
  {
      double const * const __restrict__ Qout = dhem->XoscArr[i];
      double const * const __restrict__ Pout = &dhem->XoscArr[i][3*sim->n];
      double * const __restrict__ Q_dot_out = dhem->Xosc_dotArr[i];
      double * const __restrict__ P_dot_out = &dhem->Xosc_dotArr[i][3*sim->n];

      dhem_CalculateOsculatingOrbitDerivatives_Momenta(Qout, Pout, Q_dot_out, P_dot_out);
  }
}

void dhem_Init(SIMULATION * z_sim, double z_rectificationPeriodDefault, uint32_t z_stagesPerStep)
{
  sim = z_sim;
  // Get memory for the dhem state vectors.
  dhem = (DHEM*)malloc(sizeof(DHEM));
  memset(dhem, 0, sizeof(DHEM));

  // Configure function pointers for other modules.
  sim->rhs = dhem;
  sim->f_rhs = dhem_rhs_wrapped;
  sim->fStartOfStep = dhem_CalcOscOrbitsForAllStages;
  sim->fRectify = dhem_RectifyOrbits;
  sim->fPerformSummation = dhem_PerformSummation;
  sim->fCalculateInvariant = dhem_CalculateHamiltonian;

  dhem->final_stage_index = 8;

  dhem->X = (double*)malloc(sim->stateVectorSize);
  dhem->X_dot = (double*)malloc(sim->stateVectorSize);
  dhem->rectifyTimeArray = (double*)malloc(sim->controlVectorSize);
  dhem->rectificationPeriod = (double*)malloc(sim->controlVectorSize);

  // Create space to allow for all of the osculating orbits for a step to be stored.
  dhem->XoscStore = (double*)malloc(z_stagesPerStep*sim->stateVectorSize);
  dhem->XoscArr = (double **)malloc(z_stagesPerStep*sizeof(double*));
  dhem->Xosc_dotStore = (double*)malloc(z_stagesPerStep*sim->stateVectorSize);
  dhem->Xosc_dotArr = (double **)malloc(z_stagesPerStep*sizeof(double*));
  dhem->Vosc = (double*)malloc(sim->stateVectorSize / 2);

  // Create space to allow for all of the osculating orbits for a step to be stored.
  dhem->XoscPredStore = (double*)malloc(z_stagesPerStep*sim->stateVectorSize);
  dhem->XoscPredArr = (double **)malloc(z_stagesPerStep*sizeof(double*));


  // Creat space for osculating orbit compensated summation variables
  dhem->XoscStore_cs = (double*)malloc(z_stagesPerStep*sim->stateVectorSize);
  dhem->XoscArr_cs = (double **)malloc(z_stagesPerStep*sizeof(double*));

  // Create space to allow for all of the osculating orbits for a step to be stored but in extended precision.
  dhem->XoscStore_ld = (long double*)malloc(z_stagesPerStep*sim->stateVectorSize_ld);
  dhem->XoscArr_ld = (long double **)malloc(z_stagesPerStep*sizeof(long double*));
  dhem->Xosc_dotStore_ld = (long double*)malloc(z_stagesPerStep*sim->stateVectorSize_ld);
  dhem->Xosc_dotArr_ld = (long double **)malloc(z_stagesPerStep*sizeof(long double*));
  
  // Compensated summation rhs variables
  dhem->dQdot_cs = (double*)malloc((int)sim->stateVectorSize/2);
  dhem->dQddot_cs = (double*)malloc((int)sim->stateVectorSize/2);
  dhem->dPdot_cs = (double*)malloc((int)sim->stateVectorSize/2);
  memset(dhem->dQdot_cs, 0, sim->stateVectorSize / 2);
  memset(dhem->dQddot_cs, 0, sim->stateVectorSize / 2);
  memset(dhem->dPdot_cs, 0, sim->stateVectorSize / 2);


  // Set required arrays to zero.
  memset(dhem->X, 0, sim->stateVectorSize);
  memset(dhem->X_dot, 0, sim->stateVectorSize);
  memset(dhem->rectifyTimeArray, 0, sim->controlVectorSize);
  memset(dhem->rectificationPeriod, 0, sim->controlVectorSize);

  memset(dhem->XoscStore, 0, z_stagesPerStep*sim->stateVectorSize);
  memset(dhem->XoscArr, 0, z_stagesPerStep*sizeof(double *));
  memset(dhem->Xosc_dotStore, 0, z_stagesPerStep*sim->stateVectorSize);
  memset(dhem->Xosc_dotArr, 0, z_stagesPerStep*sizeof(double *));
  memset(dhem->Vosc, 0, sim->stateVectorSize / 2);

  memset(dhem->XoscPredStore, 0, z_stagesPerStep*sim->stateVectorSize);
  memset(dhem->XoscPredArr, 0, z_stagesPerStep*sizeof(double *));

  memset(dhem->XoscStore_cs, 0, z_stagesPerStep*sim->stateVectorSize);
  memset(dhem->XoscArr_cs, 0, z_stagesPerStep*sizeof(double *));

  // Long double version.
  memset(dhem->XoscStore_ld, 0, z_stagesPerStep*sim->stateVectorSize_ld);
  memset(dhem->XoscArr_ld, 0, z_stagesPerStep*sizeof(long double *));
  memset(dhem->Xosc_dotStore_ld, 0, z_stagesPerStep*sim->stateVectorSize_ld);
  memset(dhem->Xosc_dotArr_ld, 0, z_stagesPerStep*sizeof(long double *));

  // To enable easier access to the osculating orbits.
  for(uint32_t i = 0; i < z_stagesPerStep; i++)
  {
    dhem->XoscArr[i] = &dhem->XoscStore[i*sim->stateVectorLength];
    dhem->XoscPredArr[i] = &dhem->XoscPredStore[i*sim->stateVectorLength];
    dhem->Xosc_dotArr[i] = &dhem->Xosc_dotStore[i*sim->stateVectorLength];
    dhem->XoscArr_ld[i] = &dhem->XoscStore_ld[i*sim->stateVectorLength];
    dhem->Xosc_dotArr_ld[i] = &dhem->Xosc_dotStore_ld[i*sim->stateVectorLength];    

    dhem->XoscArr_cs[i] = &dhem->XoscStore_cs[i*sim->stateVectorLength];
  }

  // Setup pointers for more human readable access.
  dhem->Qosc = dhem->XoscArr[0];
  dhem->Posc = &dhem->Qosc[3*sim->n];

  dhem->Qosc_cs = dhem->XoscArr_cs[0];
  dhem->Posc_cs = &dhem->Qosc_cs[3*sim->n];

  dhem->Qosc_dot = dhem->Xosc_dotArr[0];
  dhem->Posc_dot = &dhem->Qosc_dot[3*sim->n];

  dhem->Q = dhem->X;
  dhem->P = &dhem->X[3*sim->n];

  dhem->Q_dot = dhem->X_dot;
  dhem->P_dot = &dhem->X_dot[3*sim->n];

  dhem->m = sim->mass;
  dhem->n = sim->n;
  dhem->mTotal = 0;

  dhem->m_inv = (double*)malloc(sim->n*sizeof(double));

  for(uint32_t i = 0; i < dhem->n; i++)
  {
    dhem->m_inv[i] = 1.0 / dhem->m[i];
    dhem->mTotal += dhem->m[i];
    dhem->rectifyTimeArray[i] = sim->t0 + z_rectificationPeriodDefault;
    dhem->rectificationPeriod[i] = z_rectificationPeriodDefault;
  }

  dhem_ConvertToCOM(sim->Q0, sim->V0, sim->r0, sim->v0);
  dhem_ConvertToDHCoords(sim->Q0, sim->V0, sim->Q_dh, sim->P_dh);

  dhem_InitialiseOsculatingOrbits(sim->Q_dh, sim->P_dh, sim->t0);
}

void dhem_Free(void)
{
  // @todo clean up properly - should run a memcheck on the entire program.
  // free(dhem->X);
  // free(dhem->X_dot);
  // free(dhem->rectifyTimeArray);
  // free(dhem->rectificationPeriod);
  // free(dhem->XoscStore);
  // free(dhem->XoscArr);
  // free(dhem->Xosc_dotStore);
  // free(dhem->Xosc_dotArr);
  // free(dhem->Vosc);
  // free(dhem);
  // free(dhem->fStore);
  // free(dhem->fAccessArray);
  // free(dhem->gStore);
  // free(dhem->gAccessArray);
  // sim->rhs = NULL;
}

static inline void add_cs(double* out, double* cs, double inp)
{
    const double y = inp - cs[0];
    const double t = out[0] + y;
    cs[0] = (t - out[0]) - y;
    out[0] = t;
}
