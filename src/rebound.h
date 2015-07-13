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
#ifndef M_PI
// Make sure M_PI is defined. 
#define M_PI           3.14159265358979323846
#endif
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "integrator_sei.h"
#include "integrator_wh.h"


/**
 * Generic 3d vector
 */
struct reb_vec3d {
	double x;
	double y;
	double z;
};


/**
 * reb_particle structure.
 * @details This structure is used to represent one particle. Additional particle
 * properties should be added here. Note that when the data structure is changed, 
 * one must also update the equivalent declaration for MPI in communications_mpi.c.
 */
struct reb_particle {
	double x;	/**< x-position of the particle. */
	double y;	/**< y-position of the particle. */
	double z;	/**< z-position of the particle. */
	double vx;	/**< x-velocity of the particle. */
	double vy;	/**< y-velocity of the particle. */
	double vz;	/**< z-velocity of the particle. */
	double ax;	/**< x-acceleration of the particle. */
	double ay;	/**< y-acceleration of the particle. */
	double az;	/**< z-acceleration of the particle. */
	double m;	/**< Mass of the particle. */
#ifdef PARTICLEIDS
	int ID;		/**< Unique ID to identify particle. */
#endif
	double r; 	/**< Radius of the particle. */
	double lastcollision;	/**< Last time the particle had a physical collision. */
	struct cell* c;		/**< Pointer to the cell the particle is currently in. */
};



struct reb_context {
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
	double  max_radius[2];	/**< Two largest particle radii. Needed for collision search. */
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

	unsigned int force_is_velocitydependent; 	/**< Set to 1 if integrator needs to consider velocity dependent forces. */ 
	unsigned int gravity_ignore_10;			/**< Ignore the gravity form the central object (for WH-type integrators)*/
	struct reb_vec3d* cs;				/**< Vector containing the information for compensated gravity summation */
	int N_cs;					/**< Current number of allocated space for cs array*/


	//////////////////////////////////////////////
	/// Collisions
	struct	collision* collisions;			/**< Array of all collisions. */
	int 	collisions_NMAX;			/**< Size allocated for collisions.*/
	double coefficient_of_restitution;		/**< Used for hard sphere collision model. */
	double minimum_collision_velocity;		/**< Used for hard sphere collision model. */
	double (*collisions_coefficient_of_restitution_for_velocity) (const struct reb_context* const r, double v); /**< Used for hard sphere collision model. */
	void (*collision_resolve) (struct reb_context* const r, struct collision);
	double 	collisions_plog;			/**< Keep track of momentum exchange (used to calculate collisional viscosity in ring systems. */
	long	collisions_Nlog;			/**< Keep track of Number of collisions. */

	//////////////////////////////////////////////
	/// Variational reb_particles
	double megno_Ys;
	double megno_Yss;
	double megno_cov_Yt;	// covariance of <Y> and t
	double megno_var_t;  	// variance of t 
	double megno_mean_t; 	// mean of t
	double megno_mean_Y; 	// mean of Y
	double megno_delta0; 	// initial scale of delta (for one particle)
	long   megno_n; 	// number of covariance updates

	/**
	 * Available collision routines.
	 */
	enum {
		RB_CT_NONE = 0,
		RB_CT_DIRECT = 1,
		RB_CT_TREE = 2,
		} collision;
	/**
	 * Available integrators.
	 */
	enum {
		RB_IT_IAS15 = 0,
		RB_IT_WHFAST = 1,
		RB_IT_SEI = 2,
		RB_IT_WH = 3,
		RB_IT_LEAPFROG = 4,
		RB_IT_HYBRID = 5,
		RB_IT_NONE = 6,
		} integrator;

	/**
	 * Available boundary conditions.
	 */
	enum {
		RB_BT_NONE = 0,
		RB_BT_OPEN = 1,
		RB_BT_PERIODIC = 2,
		RB_BT_SHEAR = 3,
		} boundary;

	/**
	 * Available gravity routines.
	 */
	enum {
		RB_GT_NONE = 0,
		RB_GT_BASIC = 1,
		RB_GT_COMPENSATED = 2,
		RB_GT_TREE = 3,
		} gravity;



	//////////////////////////////////////////////
	/// reb_particles
	struct reb_particle* particles; 			/**< Main particle array. This contains all particles on this node.  */
	
	//////////////////////////////////////////////
	/// Tree
	struct cell** tree_root; 			/**< Pointer to the roots of the trees. */
	int N_tree_fixed; 				/**< reb_particle between 0 and N_tree_fixed will not be shuffled around during tree-reconstruction.  */
	double opening_angle2;	 			/**< Square of the cell opening angle \f$ \theta \f$. */

	
	//////////////////////////////////////////////
	/// Integrators
	struct reb_contextIntegratorWHFast ri_whfast;	/**< The WHFast struct */
	struct reb_contextIntegratorIAS15 ri_ias15;		/**< The IAS15 struct */
	struct reb_contextIntegratorSEI ri_sei;		/**< The SEI struct */
	struct reb_contextIntegratorWH ri_wh;		/**< The WH struct */


	//////////////////////////////////////////////
	/// Callback function
	/*
	 * This function allows the user to add additional (non-gravitational) forces.
	 */
	void (*additional_forces) (struct reb_context* const r);
	/*
	 * This function allows the user to modify the dditional (non-gravitational) forces.
	 */
	void (*post_timestep_modifications) (struct reb_context* const r);
	/**
	 * This function is called at the beginning of the simulation and at the end of
	 * each timestep.
	 */
	void (*heartbeat) (struct reb_context* r);
	/**
	 * This function is called at the end of the simulation when t>=tmax.
	 * Note that it is not called when the simulation stopped for another 
	 * reason (e.g. user interaction or crash). 
	 */ 
	void (*finished) (struct reb_context* r);

};


/**
 * Initializes all REBOUND variables and returns a REBOUND handle.. 
 * This function must be called from problem_init() before any particles are added.
 */
struct reb_context* rebound_init();

/**
 * Performon integration step.
 */
void rebound_step(struct reb_context* const r);

/**
 * Performon an integration. Starting at the current time t and until time tmax.
 * tmax==0 means integrate forever.
 */
int rebound_integrate(struct reb_context* const r, double tmax);

/**
 * Helper function to configure box.
 */
void rebound_configure_box(struct reb_context* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

/**
 * This function is called once before the integration and then after every timestep.
 * The simulation exits immediately if it returns 1.
 */
int rebound_check_exit(struct reb_context* const r, const double tmax);

/*
 * Function used to allow binary input.
 */
void rebound_reset_temporary_pointers(struct reb_context* const r);
void rebound_reset_function_pointers(struct reb_context* const r);
#endif
