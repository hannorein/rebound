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
 * This function sets up a Plummer sphere.
 * @details This function is based on a routine from the NEMO package, P. Teuben (1995).
 * For details on the implementation see the Appendix of Aarseth, Henon and Wielen (1974). 
 * @param _N Number of particles in the plummer sphere.
 * @param mlow Lower mass fraction cutoff (can be 0).
 * @param rfrac Upper radius cutoff (the Plummer sphere is formally an inifitely large object).
 * @param quiet Noisyness of the model, 0=noise, 1=medium, 2=quiet.
 * @param scale Scales the final model before adding it to the simulation.
 * @param shift Shift the entire sphere in position and velocity space (6 values). 
 */

void tools_init_plummer(int _N, double mlow, double rfrac, int quiet, double scale, double* shift);

/**
 * Initialize a particle on an orbit in the xy plane.
 * @param M Mass of the central object.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 */
struct particle tools_init_orbit2d(double M, double m, double a, double e, double omega, double f);

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
void tools_move_to_center_of_momentum();

/**
 * Returns the center of mass of particle p1 and p2.
 */
struct particle tools_get_center_of_mass(struct particle p1, struct particle p2);

#endif 	// TOOLS_H
