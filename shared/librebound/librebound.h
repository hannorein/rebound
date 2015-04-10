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

extern double 	timing;		/**< Time for last step/integration in s. */

/*
 * This functions sets the current integrator.
 * Default is IAS15. See integrator.h for options.
 */
void		integrator_set(int i);


/*
 * Integrate for one step.
 */
void rebound_step();

/* Integrate until t=_tmax.
 * The integration finisheds exactly at _tmax if 
 * exactFinishTime=1, otherwise REBOUND will overshoot slightly
 * depending on the current timestep.
 */
void integrate(double _tmax, int exactFinishTime);

#endif
