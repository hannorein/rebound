/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate particles on a circular orbit in a fixed 
 * potential.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "main.h"
#include "output.h"
#include "particle.h"
#include "boundaries.h"

void problem_init(int argc, char* argv[]){
  // Setup constants
  dt 		= M_PI*1e-3*1.1234125235345;
  boxsize 	= 2.2;
  N_active      = 2;
  init_box();

  // Initial conditions
  // The WH integrator assumes a heliocentric coordinate system.
  // Therefore the star has to be at the origin.
  struct particle star;
  star.x  = 0; star.y  = 0; star.z  = 0;
  star.vx = 0; star.vy = 0; star.vz = 0;
  star.ax = 0; star.ay = 0; star.az = 0;
  star.m  = 1;
  particles_add(star);

  // Planet 1 at 1 AU
  struct particle planet;
  planet.x  = 1.; planet.y  = 0.; planet.z  = 0.;
  double r_init1 = sqrt(planet.x*planet.x + planet.y*planet.y + planet.z*planet.z);
  planet.vx = 0;
  planet.vy = sqrt(1./r_init1);
  planet.vz = 0;
  /**/
  planet.ax = 0; planet.ay = 0; planet.az = 0;
  planet.m  = 1.e-3;
  particles_add(planet);

  // Ring particles
  double inner_R =  2./800.;
  double outer_R =  6./800.;
  int NumParticles = 20000;
  for (int ii = 1; ii <= NumParticles; ++ii) {
    double dr = 0.;
    if (NumParticles == 1)
      dr = (inner_R + outer_R)/2.;
    else
      dr = inner_R + ((ii-1.)/(NumParticles-1)) * (outer_R-inner_R);
    printf("NumParticles: %d; dr: %1.5e\n",NumParticles,dr);
    double vr = sqrt(planet.m / dr);
    struct particle ringparticle;
    ringparticle.x  = planet.x + dr;
    ringparticle.y  = 0;
    ringparticle.z  = 0;
    double r_init2 = sqrt(ringparticle.x*ringparticle.x + \
			  ringparticle.y*ringparticle.y + \
			  ringparticle.z*ringparticle.z);
    ringparticle.vx = 0;
    ringparticle.vy = sqrt(1./r_init2);
    ringparticle.vz = 0;
    ringparticle.vy = ringparticle.vy + vr;
    ringparticle.ax = 0;
    ringparticle.ay = 0;
    ringparticle.az = 0;
    ringparticle.m  = 0;
    particles_add(ringparticle);
  }
}

void problem_inloop(){
  if(output_check(1000.*dt)){
    output_timing();
  }
  if(output_check(124.)){
    //output_append_ascii("position.txt");
    output_append_orbits("orbits.txt");
  }
}

void problem_output(){
}

void problem_finish(){
}
