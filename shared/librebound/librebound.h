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
extern int 	root_nx;	/**< Number of root boxes in x direction. Default: 1. */
extern int 	root_ny;	/**< Number of root boxes in y direction. Default: 1. */
extern int 	root_nz;	/**< Number of root boxes in z direction. Default: 1. */
extern double 	boxsize;	/**< Size of a root box. Needs to be set in problem_init(). */
extern double 	boxsize_x;	/**< Size of the entire box in the x direction, root_nx*boxsize. Set in init_box().*/
extern double 	boxsize_y;	/**< Size of the entire box in the y direction, root_ny*boxsize. Set in init_box().*/
extern double 	boxsize_z;	/**< Size of the entire box in the z direction, root_nz*boxsize. Set in init_box().*/

extern double 	timing_initial;	/**< System time at start. Used to meassure total cpu time. */

// Integer flag that determines what kind of integrator is used.
extern int selected_integrator; 

#endif
