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
#include <math.h>
#include "Simulation.h"
#include "dhem.h"
#include "radau.h"
#include "radau_step.h"
#include "UniversalVars.h"
#include "terminate.h"
#include <stdlib.h>
#include <string.h>

static uint32_t Intersection(uint32_t * bodies);
static uint32_t Collision(void);
static double ClosestApproach(uint32_t * z_bodies);
static uint32_t CombinedDataRecording(void);
static uint32_t IntersectionWrapped(void);

SIMULATION * sim = NULL;

double * Q = NULL;
double * P = NULL;

double cutOffForRecording = 1E-2;
static uint32_t FINAL_STAGE_INDEX = 8;

void Terminate_Init(SIMULATION * z_sim, uint32_t (**f)(void), TERMINATION terminate)
{
  sim = z_sim;

  if(terminate == none)
  {
    f = NULL;
  }
  else if (terminate == intersection) {
    *f = IntersectionWrapped;
  }
  else if (terminate == collision)
  {
    *f = Collision;
    Q = (double*)malloc(3*sim->n*sizeof(double));
    P = (double*)malloc(3*sim->n*sizeof(double));
  }
  else if(terminate == combined)
  {
    *f = CombinedDataRecording;
    Q = (double*)malloc(3*sim->n*sizeof(double));
    P = (double*)malloc(3*sim->n*sizeof(double));
  }

  char * encounter_file_string = (char *)malloc(strlen(sim->outputFile)+strlen("_encounter.txt")+1);
  strcpy(encounter_file_string, sim->outputFile);
  uint32_t strLen = strlen(encounter_file_string);
  encounter_file_string[strLen-4] = 0;      // Remove .txt from output file string [probs not the best way of doing this incase the extension changes.]
  strcat(encounter_file_string, "_encounter.txt");
  
  sim->encounterFile = fopen(encounter_file_string, "w+");
  free(encounter_file_string);  
}


void Terminate_Free(SIMULATION * z_sim)
{
  free(Q);
  free(P);
}

static uint32_t IntersectionWrapped(void)
{
  uint32_t bodies[2];
  return Intersection(bodies);
}

static uint32_t Intersection(uint32_t * bodies)
{
  uint32_t retVal = 0;

  for(uint32_t i = 1; i < sim->n-1; i++)
  {
    if(sim->uVars->apo[i] > sim->uVars->peri[i+1])
    {
      retVal = 1;
      bodies[0] = i;
      bodies[1] = i+1;
      break;
    }
  }
  return retVal;
}

double cutoff = (2*6378)/1.495978707E8;

static uint32_t Collision(void)
{
    uint32_t retVal = 0;
    double approach = 0;
    uint32_t bodies[2] = {0,0};

    approach = ClosestApproach(bodies);

    if(approach < cutoff)
    {
      retVal = 1;
    }

  return retVal;
}

static double ClosestApproach(uint32_t * z_bodies)
{
  double sep[3] = {0,0,0};
  double sep_norm = 0.0;
  double closest = 1E10;

  sim->fPerformSummation(Q, P, sim->radau->dQ, sim->radau->dP, FINAL_STAGE_INDEX);

  double * r = Q;

  for(uint32_t i = 0; i < sim->n; i++)
  {
    for(uint32_t j = i+1; j < sim->n; j++)
    {
      sep[0] = r[3*i] - r[3*j];
      sep[1] = r[3*i+1] - r[3*j+1];
      sep[2] = r[3*i+2] - r[3*j+2];
      sep_norm = sqrt(sep[0]*sep[0]+sep[1]*sep[1]+sep[2]*sep[2]);

      if(sep_norm < closest)
      {
        closest = sep_norm;
        z_bodies[0] = i;
        z_bodies[1] = j;
      }
    }
  }
  return closest;
}

static uint32_t CombinedDataRecording(void)
{
  static uint32_t intersectionOccurred = 0;
  static double closest_approach = 1E30;
  static double approach_last = 1E30;
  static uint32_t output_flag = 0;
  double approach = 0;
  uint32_t retVal = 0;
  uint32_t bodies[2] = {0,0};

  if(intersectionOccurred == 0)
  {
    intersectionOccurred = Intersection(bodies);

    if(intersectionOccurred != 0)
    {
      sim->fPerformSummation(sim->radau->Qout, sim->radau->Pout, sim->radau->dQ, sim->radau->dP, FINAL_STAGE_INDEX);
      double H = sim->fCalculateInvariant(sim->radau->Qout, sim->radau->Pout);
      H = fabs((H-sim->H0)/sim->H0);
      fprintf(sim->encounterFile, "%.16E %.16E 0.0 0.0 0.0 0.0 0.0 0.0 0.0 %d %d\n", sim->radau->t, H, bodies[0], bodies[1]);
    }
  }

  approach = ClosestApproach(bodies);

  if(approach < closest_approach)
  {
      closest_approach = approach;
      output_flag = 1;
  }

  // Dont record the closest point until we start moving away (otherwise we get spammed on the way in)
  if((approach > approach_last && output_flag) || closest_approach < cutoff)
  {      
      output_flag = 0;
    
      sim->fPerformSummation(sim->radau->Qout, sim->radau->Pout, sim->radau->dQ, sim->radau->dP, FINAL_STAGE_INDEX);
      double H = sim->fCalculateInvariant(sim->radau->Qout, sim->radau->Pout);
      H = fabs((H-sim->H0)/sim->H0);     
    
      double xi, yi, vxi, vyi, vxj, vyj;
      uint32_t i = bodies[0];
      uint32_t j = bodies[1];

      xi = sim->radau->Qout[3*i+0];
      yi = sim->radau->Qout[3*i+1];

      vxi = sim->radau->Pout[3*i+0]/sim->mass[i];
      vyi = sim->radau->Pout[3*i+1]/sim->mass[i];

      vxj = sim->radau->Pout[3*j+0]/sim->mass[j];
      vyj = sim->radau->Pout[3*j+1]/sim->mass[j];

      fprintf(sim->encounterFile, "%.16E %.16E %.16E %.3E %.3E %.3E %.3E %.3E %.3E %d %d\n", 
              sim->radau->t, H, closest_approach, xi, yi, vxi, vyi, vxj, vyj, bodies[0], bodies[1]);   

      // This is the distance we class as a collision.
      if(closest_approach < cutoff)
      {
        retVal = 1;
      }                  
  }
  approach_last = approach;

  return retVal;
}
