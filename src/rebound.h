/**
 * @file 	main.h
 * @brief 	Main header file with widely used global variables.
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
#ifndef _MAIN_H
#define _MAIN_H
#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)
#define TREE
#endif
#ifndef M_PI
// Make sure M_PI is defined. 
#define M_PI           3.14159265358979323846
#endif
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "integrator_sei.h"
#include "integrator_wh.h"

/*
 * Available integrators.
 */
typedef enum {
	IAS15 = 0,
	WHFAST = 1,
	SEI = 2,
	WH = 3,
	LEAPFROG = 4,
	HYBRID = 5,
	NONE = 6,
	} integrator_t;

struct Rebound {
	double 	t;		/**< Current simulation time. */
	double 	G;		/**< Gravitational constant. Default: 1. */
	double 	softening;	/**< Gravitational softening parameter. Default: 0. */
	double 	dt;		/**< Current timestep. */
	double 	dt_last_done;	/**< Last full timestep (used if exact_finish_time==1). */
	double 	boxsize;	/**< Size of a root box. Needs to be set in problem_init(). */
	double 	boxsize_x;	/**< Size of the entire box in the x direction, root_nx*boxsize. Set in box_init().*/
	double 	boxsize_y;	/**< Size of the entire box in the y direction, root_ny*boxsize. Set in box_init().*/
	double 	boxsize_z;	/**< Size of the entire box in the z direction, root_nz*boxsize. Set in box_init().*/
	double 	boxsize_max;	/**< Maximum size of the entire box in any direction. Set in box_init().*/
	int 	N;		/**< Current number of particles on this node. */
	int 	Nmax;		/**< Current maximum space allocated in the particles array on this node. */
	int 	N_active;	/**< Number of massive particles included in force calculation. Default: N.*/
	int 	N_megno;	/**< Number of megno particles. Default: 0.*/
	int 	root_nx;	/**< Number of root boxes in x direction. Default: 1. */
	int 	root_ny;	/**< Number of root boxes in y direction. Default: 1. */
	int 	root_nz;	/**< Number of root boxes in z direction. Default: 1. */
	int 	root_n;		/**< Total number of root boxes in all directions, root_nx*root_ny*root_nz. Default: 1. Set in box_init().*/
	int	nghostx;	/**< Number of ghostboxes in x direction. */
	int 	nghosty;	/**< Number of ghostboxes in y direction. */
	int 	nghostz;	/**< Number of ghostboxes in z direction. */
	int 	exit_simulation;/**< Set to 1 to exit the simulation at the end of the next timestep. */
	int 	exact_finish_time; /**< Set to 1 to finish the integration exactly at tmax. Set to 0 to finish at the next dt. */
	integrator_t integrator; /**< Variable setting the current integrator.  */

	unsigned int force_is_velocitydependent; 	/**< Set to 1 if integrator needs to consider velocity dependent forces. */ 
	unsigned int gravity_ignore_10;			/**< Ignore the gravity form the central object (for WH-type integrators)*/

	//////////////////////////////////////////////
	/// Variational Particles
	double megno_Ys;
	double megno_Yss;
	double megno_cov_Yt;	// covariance of <Y> and t
	double megno_var_t;  	// variance of t 
	double megno_mean_t; 	// mean of t
	double megno_mean_Y; 	// mean of Y
	double megno_delta0; 	// initial scale of delta (for one particle)
	long   megno_n; 	// number of covariance updates


	//////////////////////////////////////////////
	/// Particles
	struct Particle* particles; 			/**< Main particle array. This contains all particles on this node.  */
	
#ifdef TREE
	//////////////////////////////////////////////
	/// Tree
	struct cell** tree_root; 			/**< Pointer to the roots of the trees. */
	int N_tree_fixed; 				/**< Particle between 0 and N_tree_fixed will not be shuffled around during tree-reconstruction.  */
	double opening_angle2;	 			/**< Square of the cell opening angle \f$ \theta \f$. */

#endif // TREE
	
	//////////////////////////////////////////////
	/// Integrators
	struct ReboundIntegratorWHFast ri_whfast;	/**< The WHFast struct */
	struct ReboundIntegratorIAS15 ri_ias15;		/**< The IAS15 struct */
	struct ReboundIntegratorSEI ri_sei;		/**< The SEI struct */
	struct ReboundIntegratorWH ri_wh;		/**< The WH struct */


	//////////////////////////////////////////////
	/// Callback function
	/*
	 * This function allows the user to add additional (non-gravitational) forces.
	 */
	void (*additional_forces) (struct Rebound* const r);
	/*
	 * This function allows the user to modify the dditional (non-gravitational) forces.
	 */
	void (*post_timestep_modifications) (struct Rebound* const r);
	/**
	 * This function is called at the beginning of the simulation, at the end of
	 * each timestep and at the end of the simulation. 
	 */
	void (*post_timestep) (struct Rebound* r);
	/**
	 * This function is called at the end of the simulation when t>=tmax.
	 * Note that it is not called when the simulation stopped for another 
	 * reason (e.g. user interaction or crash). 
	 */ 
	void (*finished) (struct Rebound* r);

};


/**
 * Initializes all REBOUND variables and returns a REBOUND handle.. 
 * This function must be called from problem_init() before any particles are added.
 */
struct Rebound* rebound_init();

/**
 * Performon integration step.
 */
void rebound_step(struct Rebound* const r);

/**
 * Performon an integration. Starting at the current time t and until time tmax.
 */
int rebound_integrate(struct Rebound* const r, double tmax);

/**
 * Helper function to configure box.
 */
void rebound_configure_box(struct Rebound* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

/**
 * This function is called once before the integration and then after every timestep.
 * The simulation exits immediately if it returns 1.
 */
int rebound_check_exit(struct Rebound* const r, const double tmax);
#endif
