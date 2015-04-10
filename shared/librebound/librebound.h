#ifndef _LIBREBOUND_H
#define _LIBREBOUND_H

#ifndef M_PI
#define M_PI 3.1415926535879323846
#endif

extern double 	softening;	/**< Gravitational softening parameter. Default: 0. */
extern double 	G;		/**< Gravitational constant. Default: 1. */
extern double 	t;		/**< Current simulation time. */
extern double 	tmax;		/**< Maximum simulation time. Simulation stops if t>=tmax. Simulation runs forever if t==0.*/
extern double 	dt;		/**< Current timestep. */
extern int 	exit_simulation;/**< Set to 1 to exit the simulation at the end of the timestep. */
extern int 	N;		/**< Current number of particles on this node. Ns set in particle.c*/
extern int 	N_active;	/**< Number of massive particles included in force calculation. Default: N.*/
extern int 	N_megno;	/**< Number of megno particles. Default: 0.*/

extern double 	timing_initial;	/**< System time at start. Used to meassure total cpu time. */

/*
 * This functions sets the current integrator.
 * Default is IAS15. See integrator.h for options.
 */
void		integrator_set(int i);

#endif
