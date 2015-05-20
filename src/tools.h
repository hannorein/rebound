/**
 * @file 	tools.h
 * @brief 	Tools for creating distributions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#ifndef TOOLS_H
#define TOOLS_H
#include "particle.h"
/**
 * Struct representing a Keplerian orbit.
 */
struct orbit {
	double a;
	double r;	// Radial distance from central object
	double h;	// Angular momentum
	double P;	// Orbital period
	double l;
	double e;
	double inc;
	double Omega; 	// longitude of ascending node
	double omega; 	// argument of perihelion
	double f; 	// true anomaly
};


/**
 * Calculates a random variable in a given range.
 * @param min Minimum value.
 * @param max Maximum value.
 */
double tools_uniform(double min, double max);

/**
 * Calculates a random variable drawn form a powerlaw distribution.
 * @param min Minimum value.
 * @param max Maximum value.
 * @param slop Slope of powerlaw distribution.
 */
double tools_powerlaw(double min, double max, double slope);

/**
 * Calculate a random number with normal distribution.
 * Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
 * @param variance Variance of normal distribution.
 * @return Random number with normal distribution (mean 0). 
 */
double tools_normal(double variance);

/**
 * Calculates a random variable drawn form a Rayleigh distribution.  Calculated as described on Rayleigh distribution wikipedia page
 * @param sigma Scale parameter.
 */
double tools_rayleigh(double sigma);

/**
 * This function sets up a Plummer sphere.
 * @param _N Number of particles in the plummer sphere.
 * @param M Total mass of the cluster.
 * @param R Characteristic radius of the cluster.
 */

void tools_init_plummer(int _N, double M, double R);

/**
 * Initialize a particle on an orbit in the xy plane.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 * @param f true anomaly of the particle.
 */
struct particle tools_init_orbit2d(double M, double m, double a, double e, double omega, double f);

/**
 * Initialize a particle on a 3D orbit.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane.
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 */

struct particle tools_init_orbit3d(double M, double m, double a, double e, double i, double Omega, double omega, double f);

/**
 * This function calculated orbital elements for a given particle. 
 * @param p Particle for which the orbit is calculated.
 * @param star Star or central object particle
 * @return Orbital parameters. 
 */
struct orbit tools_p2orbit(struct particle p, struct particle star);

/**
 * Move to center of momentum and center of mass frame.
 */
void tools_move_to_center_of_momentum(void);

/**
 * Returns the center of mass of particle p1 and p2.
 */
struct particle tools_get_center_of_mass(struct particle p1, struct particle p2);

/* 
 * Init the MEGNO particles
 **/
void tools_megno_init(double delta);

/*
 * Returns the current value of <Y>
 **/
double tools_megno(void);

/*
 * Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
 **/
double tools_lyapunov(void);

/*
 * Returns deltad/delta (Note, there is a typo in Gozdziewski et al 2001).
 **/

double tools_megno_deltad_delta(void);

/*
 * Update MEGNO after a successful timestep by adding dY (=ddelta/delta*dt)
 **/
void tools_megno_update(double dY);

/**
 * Calculate the total energy (potential and kinetic).
 * Might not work for WH.
 * @return Total energy. 
 */
double tools_energy(void);

#endif 	// TOOLS_H
