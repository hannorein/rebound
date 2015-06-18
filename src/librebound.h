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

extern const char *build_str;	/**< Contains last compile date/time information. */

extern int closeEncounterPi; 	/**< IDs of the particles which had a close encounter */
extern int closeEncounterPj;

/*
 * This functions sets the current integrator.
 * Default is IAS15. See integrator.h for options.
 */
void		integrator_set(int i);


/*
 * Integrate for one step.
 */
void rebound_step(int do_timing);

/* Integrate until t=_tmax.
 * The integration finisheds exactly at _tmax if 
 * exactFinishTime=1, otherwise REBOUND will overshoot slightly
 * depending on the current timestep.
 * If the integrator_whfast_safe_mode flag is set to 1 (default), then the integrator
 * will synchronize the positions and velocities after every timestep.
 * This will cause the integrator to be slower and less accurate. To turn off 
 * safe_mode, and take the appropriate steps manually, see AdvWHFast.ipynb
 * in the python_tutorials folder (navigate to it on github
 * if you don't have ipython notebook installed).  The explanation is general, and
 * the python and C flags have the same names.
 *
 * If maxR or minD are set to a value other than 0, then after every
 * timetep, REBOUND checks if a particle escaped or if two  particles
 * collided. Note that this can be a slowdown
 *
 * Return values:
 *   0 = All good
 *   1 = No particles left
 *   2 = Particle distance exceeds maxR
 *   3 = Close encounter closer than minD
 */
int integrate(double _tmax, int exact_finish_time, double maxR, double minD);

/*
 * This function allows the user to add additional (non-gravitational) forces.
 */
extern void (*problem_additional_forces) (void);

/*
 * This function allows the user to add additional (non-gravitational) forces.
 */
extern void (*problem_additional_forces_with_parameters) (struct particle* particles, double t, double dt, double G, int N, int N_megno);

/*
 * This function allows the user to modify particles after a timestep is completed.
 */

extern void (*post_timestep_modifications) (void);
/*
 * This function allows the user to modify particles after a timestep is completed. 
 */
extern void (*post_timestep_modifications_with_parameters) (struct particle* particles, double t, double dt, double G, int N, int N_megno);

#endif
