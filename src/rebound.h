/**
 * @file    rebound.h
 * @brief   REBOUND API definition.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
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
#define M_PI           3.14159265358979323846       ///< The mathematical constant pi.
#endif
#ifdef MPI
#include "mpi.h"
#endif // MPI

extern const char* reb_build_str;   ///< Date and time build string.
extern const char* reb_version_str; ///< Version string.

/**
 * @brief Enumeration describing the return status of rebound_integrate
 */
enum REB_STATUS {
    REB_RUNNING_PAUSED = -3,    ///< Simulation is paused by visualization.
    REB_RUNNING_LAST_STEP = -2, ///< Current timestep is the last one. Needed to ensures that t=tmax exactly.
    REB_RUNNING = -1,           ///< Simulation is current running, no error occured.
    REB_EXIT_SUCCESS = 0,       ///< Integration finished successfully.
    REB_EXIT_ERROR = 1,     ///< A generic error occured and the integration was not successfull.
    REB_EXIT_NOPARTICLES = 2,   ///< The integration ends early because no particles are left in the simulation.
    REB_EXIT_ENCOUNTER = 3,     ///< The integration ends early because two particles had a close encounter (see exit_min_distance)
    REB_EXIT_ESCAPE = 4,        ///< The integration ends early because a particle escaped (see exit_max_distance)  
    REB_EXIT_USER = 5,      ///< User caused exit, simulation did not finish successfully.
};

struct reb_simulation;

/**
 * @brief Generic 3d vector, for internal use only.
 */
struct reb_vec3d {
    double x; ///< x coordinate
    double y; ///< y coordinate
    double z; ///< z coordinate
};

/**
 * @brief Generic 7d pointer, for internal use only (IAS15).
 */
struct reb_dp7 {
    double* restrict p0; ///< 0 substep
    double* restrict p1; ///< 1 substep
    double* restrict p2; ///< 2 substep
    double* restrict p3; ///< 3 substep
    double* restrict p4; ///< 4 substep
    double* restrict p5; ///< 5 substep
    double* restrict p6; ///< 6 substep
};

/**
 * @details Structure that contains the relative position and velocity of a ghostbox.
 */
struct reb_ghostbox{
    double shiftx;      ///< Relative x position
    double shifty;      ///< Relative y position
    double shiftz;      ///< Relative z position
    double shiftvx;     ///< Relative x velocity
    double shiftvy;     ///< Relative y velocity
    double shiftvz;     ///< Relative z velocity
};


/**
 * @brief Structure representing one REBOUND particle.
 * @details This structure is used to represent one particle. 
 * If this structure is changed, the corresponding python structure
 * needs to be changes as well. Also update the equivalent declaration 
 * for MPI in communications_mpi.c.
 */
struct reb_particle {
    double x;           ///< x-position of the particle. 
    double y;           ///< y-position of the particle. 
    double z;           ///< z-position of the particle. 
    double vx;          ///< x-velocity of the particle. 
    double vy;          ///< y-velocity of the particle. 
    double vz;          ///< z-velocity of the particle. 
    double ax;          ///< x-acceleration of the particle. 
    double ay;          ///< y-acceleration of the particle. 
    double az;          ///< z-acceleration of the particle. 
    double m;           ///< Mass of the particle. 
    double r;           ///< Radius of the particle. 
    double lastcollision;       ///< Last time the particle had a physical collision.
    struct reb_treecell* c;     ///< Pointer to the cell the particle is currently in.
    int id;             ///< Unique id to identify particle.
    void* ap;           ///< Functionality for externally adding additional properties to particles.
    struct reb_simulation* sim; ///< Pointer to the parent simulation.
};


/**
 * @brief Structure representing a Keplerian orbit.
 * @details This structure is returned when calculating 
 * a Keplerian orbit from Cartesian coordinates. 
 */
struct reb_orbit {
    double d;        ///< Radial distance from central object
    double v;        ///< velocity relative to central object's velocity
    double h;        ///< Angular momentum
    double P;        ///< Orbital period
    double n;        ///< Mean motion
    double a;        ///< Semi-major axis
    double e;        ///< Eccentricity
    double inc;      ///< Inclination
    double Omega;    ///< Longitude of ascending node
    double omega;    ///< Argument of pericenter
    double pomega;   ///< Longitude of pericenter
    double f;        ///< True anomaly
    double M;        ///< Mean anomaly
    double l;        ///< Mean Longitude
    double theta;    ///< True Longitude
    double T;        ///< Time of pericenter passage
};


////////////////////////////////
// Integrator structs 

/**
 * @brief This structure contains variables and pointer used by the HYBRID integrator.
 */
struct reb_simulation_integrator_hybrid {
    double switch_ratio;    ///< Default is 8 mutual Hill Radii 
    enum {
        SYMPLECTIC,     ///< HYBRID integrator is currently using a symplectic integrator
        HIGHORDER   ///< HYBRID integrator is currently using a high order non-symplectic integrator
        } 
        mode;       ///< Flag determining the current integrator used
};

/**
 * @brief This structure contains variables and pointer used by the IAS15 integrator.
 */
struct reb_simulation_integrator_ias15 {
    /**
     * @brief This parameter controls the accuracy of the integrator.
     * @details Set to 0 to make IAS15 a non-adaptive integrator.
     * The default value is: 1e-9.
     **/
    double epsilon;

    /**
     * @brief The minimum allowed timestep.
     * @details The default value is 0 (no minimal timestep).
     * Set a finite value to this variable if the IAS15 integrator has problems
     * and the timestep becomes excessively small.
     **/
    double min_dt;
    
    /** 
     * @brief Flag that determines how relative acceleration error is estimated.
     * @details If set to 1, estimate the fractional error by max(acceleration_error)/max(acceleration), 
     * where max is take over all particles. If set to 0, estimate the fractional error by 
     * max(acceleration_error/acceleration).
     **/
    unsigned int epsilon_global;


    
    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    /**
     * @brief Counter how many times the iteration did not converge. 
     */
    unsigned long iterations_max_exceeded;



    int allocatedN;             ///< Size of allocated arrays.

    double* restrict at;            ///< Temporary buffer for acceleration
    double* restrict x0;            ///<                      position (used for initial values at h=0) 
    double* restrict v0;            ///<                      velocity
    double* restrict a0;            ///<                      acceleration
    double* restrict csx;           ///<                      compensated summation for x
    double* restrict csv;           ///<                      compensated summation for v
    double* restrict csa0;          ///<                      compensated summation for a

    struct reb_dp7 g;
    struct reb_dp7 b;
    struct reb_dp7 csb;         ///< Compensated summation for b
    struct reb_dp7 e;

    // The following values are used for resetting the b and e coefficients if a timestep gets rejected
    struct reb_dp7 br;
    struct reb_dp7 er;
    /**
     * @endcond
     */

};

/**
 * @brief This structure contains variables used by the SEI integrator.
 * @details This is where the user sets the orbital frequency OMEGA for 
 * shearing sheet simulations.
 */
struct reb_simulation_integrator_sei {
    double OMEGA;       ///< Epicyclic/orbital frequency.
    double OMEGAZ;      ///< Epicyclic frequency in vertical direction.

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    double lastdt;      ///< Cached sin(), tan() for this value of dt.
    double sindt;       ///< Cached sin() 
    double tandt;       ///< Cached tan() 
    double sindtz;      ///< Cached sin(), z axis
    double tandtz;      ///< Cached tan(), z axis
    /** @endcond */
};

/**
 * @brief This structure contains variables used by the WH integrator.
 * @details Nothing needs to be changed by the user. All the variables are just for internal use.
 */
struct reb_simulation_integrator_wh {
    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    int allocatedN;
    double* eta;
    /** @endcond */
};

/**
 * @brief This structure contains variables used by the WHFast integrator.
 */
struct reb_simulation_integrator_whfast {
    /**
     * @brief This variable turns on/off different symplectic correctors for WHFast.
     * @details 
     * - 0 (default): turns off all correctors
     * - 3: uses third order (two-stage) corrector 
     * - 5: uses fifth order (four-stage) corrector 
     * - 7: uses seventh order (six-stage) corrector 
     * - 11: uses eleventh order (ten-stage) corrector 
     */
    unsigned int corrector;

    /** 
     * @brief Setting this flag to one will recalculate Jacobi coordinates from the particle structure in the next timestep. 
     * @details After the timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_jacobi_this_timestep;

    /**
     * @brief If this flag is set (the default), whfast will recalculate jacobi coordinates and synchronize
     * every timestep, to avoid problems with outputs or particle modifications
     * between timesteps. 
     * @details Setting it to 0 will result in a speedup, but care
     * must be taken to synchronize and recalculate jacobi coordinates when needed.
     * See AdvWHFast.ipynb in the python_tutorials folder (navigate to it on github
     * if you don't have ipython notebook installed).  The explanation is general, and
     * the python and C flags have the same names.
     */
    unsigned int safe_mode;

    /**
     * @brief Jacobi coordinates
     * @details This array contains the Jacobi coordinates of all particles.
     * It is automatically filled and updated by WHfast.
     * Access this array with caution.
     */
    struct reb_particle* restrict p_j;

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    double* restrict eta;       ///< Struct containg Jacobi eta parameters 
    double Mtotal;          ///< Total mass, used for Jacobi coordinates 

    unsigned int is_synchronized;   ///< Flag to determine if current particle structure is synchronized
    unsigned int allocated_N;   ///< Space allocated in arrays
    unsigned int timestep_warning;  ///< Counter of timestep warnings
    unsigned int recalculate_jacobi_but_not_synchronized_warning;   ///< Counter of Jacobi synchronization errors
    /**
     * @endcond
     */
};


/**
 * @brief Collision structure describing a single collision.
 * @details This structure is used to save a collision during collision search. 
 * It is passed to the collision_resolve function.
 */
struct reb_collision{
    int p1;         ///< One of the colliding particles
    int p2;         ///< One of the colliding particles
    struct reb_ghostbox gb; ///< Ghostbox (of particle p1, used for periodic and shearing sheet boundary conditions)
#if defined(COLLISIONS_SWEEP) || defined(COLLISIONS_SWEEPPHI)
    double time;        ///< Time of collision.
    int crossing;       ///< Collision occurs at the interface of two sweep boxes.
#endif // COLLISIONS_SWEEP
    int ri;         ///< Index of rootcell (needed for MPI only).
};


/**
 * @brief Struct describing the properties of a set of variational equations.
 * @details One struct describes one or more sets of variational equations.
 * If testparticle is set to -1, then it is assumed that all particles are massive
 * and all particles influence all other particles. If testparticle is >=0 then 
 * the particle with that index is assumed to be a testparticle, i.e. it does not 
 * influence other particles. For second order variational equation, index_1st_order_a/b 
 * is the index in the particle array that corresponds to the 1st order variational 
 * equations.
 */
struct reb_variational_configuration{
    struct reb_simulation* sim; ///< Reference to the simulation.
    int order;                  ///< Order of the variational equation. 1 or 2. 
    int index;                  ///< Index of the first variational particle in the particles array.
    int testparticle;           ///< Is this variational configuration describe a test particle? -1 if not.
    int index_1st_order_a;      ///< Used for 2nd order variational particles only: Index of the first first order variational particle in the particles array.
    int index_1st_order_b;      ///< Used for 2nd order variational particles only: Index of the first first order variational particle in the particles array.
};


/**
 * @brief Main struct encapsulating one entire REBOUND simulation
 * @details This structure contains all variables, status flags and pointers of one 
 * REBOUND simulation. To create a REBOUND simulation use the reb_create_simulation()
 * function. This will ensure that all variables and pointers are initialized correctly.
 */
struct reb_simulation {
    /**
     * \name Variables related to time, current number of particles and simulation status/control 
     * @{
     */
    double  t;                      ///< Current simulation time. 
    double  G;                      ///< Gravitational constant. Default: 1. 
    double  softening;              ///< Gravitational softening parameter. Default: 0. 
    double  dt;                     ///< Current timestep. 
    double  dt_last_done;           ///< Last dt used by integrator
    int     N;                      ///< Current number of particles on this node. 
    int     N_var;                  ///< Total number of variational particles. Default: 0.
    int     var_config_N;           ///< Number of variational configuration structs. Default: 0.
    struct reb_variational_configuration* var_config;   ///< These configuration structs contain details on variational particles. 
    int     N_active;               ///< Number of massive particles included in force calculation (default: N). Particles with index >= N_active are considered testparticles.
    int     testparticle_type;      ///< Type of the particles with an index>=N_active. 0 means particle does not influence any other particle (default), 1 means particles with index < N_active feel testparticles (similar to MERCURY's small particles). Testparticles never feel each other.
    int     allocatedN;             ///< Current maximum space allocated in the particles array on this node. 
    struct reb_particle* particles; ///< Main particle array. This contains all particles on this node.  
    struct reb_vec3d* gravity_cs;   ///< Vector containing the information for compensated gravity summation 
    int     gravity_cs_allocatedN;  ///< Current number of allocated space for cs array
    struct reb_treecell** tree_root;///< Pointer to the roots of the trees. 
    int     tree_needs_update;      ///< Flag to force a tree update (after boundary check)
    double opening_angle2;          ///< Square of the cell opening angle \f$ \theta \f$. 
    enum REB_STATUS status;         ///< Set to 1 to exit the simulation at the end of the next timestep. 
    int     exact_finish_time;      ///< Set to 1 to finish the integration exactly at tmax. Set to 0 to finish at the next dt. Default is 1. 

    unsigned int force_is_velocity_dependent;   ///< Set to 1 if integrator needs to consider velocity dependent forces.  
    unsigned int gravity_ignore_10; ///< Ignore the gravity form the central object (for WH-type integrators)
    double output_timing_last;      ///< Time when reb_output_timing() was called the last time. 
    double exit_max_distance;       ///< Exit simulation if distance from origin larger than this value 
    double exit_min_distance;       ///< Exit simulation if distance from another particle smaller than this value 
    double usleep;                  ///< Wait this number of microseconds after each timestep, useful for slowing down visualization. Set to negative value to disable visualization (despite compiling with OPENGL=1).  
    /** @} */

    /**
     * \name Variables related to ghost/root boxes 
     * @{
     */
    struct  reb_vec3d boxsize;  ///< Size of the entire box, root_x*boxsize. 
    double  boxsize_max;        ///< Maximum size of the entire box in any direction. Set in box_init().
    double  root_size;      ///< Size of a root box. 
    int     root_n;         ///< Total number of root boxes in all directions, root_nx*root_ny*root_nz. Default: 1. Set in box_init().
    int     root_nx;        ///< Number of root boxes in x direction. Default: 1. 
    int     root_ny;        ///< Number of root boxes in y direction. Default: 1. 
    int     root_nz;        ///< Number of root boxes in z direction. Default: 1. 
    int     nghostx;        ///< Number of ghostboxes in x direction. 
    int     nghosty;        ///< Number of ghostboxes in y direction. 
    int     nghostz;        ///< Number of ghostboxes in z direction. 
    /** @} */
#ifdef MPI
    /**
     * \name Variables related to MPI 
     * @{
     */
    int    mpi_id;                              ///< Unique id of this node (starting at 0). Used for MPI only.
    int    mpi_num;                             ///< Number of MPI nodes. Used for MPI only.
    MPI_Datatype mpi_particle;                  ///< MPI datatype corresponding to the C struct reb_particle. 
    struct reb_particle** particles_send;       ///< Send buffer for particles. There is one buffer per node. 
    int*   particles_send_N;                    ///< Current length of particle send buffer. 
    int*   particles_send_Nmax;                 ///< Maximal length of particle send beffer before realloc() is needed. 
    struct reb_particle** particles_recv;       ///< Receive buffer for particles. There is one buffer per node. 
    int*   particles_recv_N;                    ///< Current length of particle receive buffer. 
    int*   particles_recv_Nmax;                 ///< Maximal length of particle receive beffer before realloc() is needed. */

    MPI_Datatype mpi_cell;                      ///< MPI datatype corresponding to the C struct reb_treecell. 
    struct reb_treecell** tree_essential_send;  ///< Send buffer for cells. There is one buffer per node. 
    int*   tree_essential_send_N;               ///< Current length of cell send buffer. 
    int*   tree_essential_send_Nmax;            ///< Maximal length of cell send beffer before realloc() is needed. 
    struct reb_treecell** tree_essential_recv;  ///< Receive buffer for cells. There is one buffer per node. 
    int*   tree_essential_recv_N;               ///< Current length of cell receive buffer. 
    int*   tree_essential_recv_Nmax;            ///< Maximal length of cell receive beffer before realloc() is needed. 
    /** @} */
#endif // MPI

    /**
     * \name Variables related to collision search and detection 
     * @{
     */
    struct reb_collision* collisions;       ///< Array of all collisions. 
    int collisions_allocatedN;          ///< Size allocated for collisions.
    double minimum_collision_velocity;      ///< Used for hard sphere collision model. 
    double collisions_plog;             ///< Keep track of momentum exchange (used to calculate collisional viscosity in ring systems. 
    double max_radius[2];               ///< Two largest particle radii, set automatically, needed for collision search. 
    long collisions_Nlog;               ///< Keep track of number of collisions. 
    /** @} */

    /**
     * \name Variables related to the chaos indicator MEGNO 
     * @{
     */
    int calculate_megno;    ///< Internal flag that determines if megno is calculated (default=0, but megno_init() sets it to the index of variational particles used for megno)
    double megno_Ys;    ///< Running megno sum (internal use)
    double megno_Yss;   ///< Running megno sum (internal use)
    double megno_cov_Yt;    ///< covariance of MEGNO Y and t
    double megno_var_t;     ///< variance of t 
    double megno_mean_t;    ///< mean of t
    double megno_mean_Y;    ///< mean of MEGNO Y
    long   megno_n;     ///< number of covariance updates
    /** @} */

    /**
     * \name Variables describing the current module selection 
     * @{
     */
    /**
     * @brief Available collision routines
     */
    enum {
        REB_COLLISION_NONE = 0,     ///< Do not search for collisions (default)
        REB_COLLISION_DIRECT = 1,   ///< Direct collision search O(N^2)
        REB_COLLISION_TREE = 2,     ///< Tree based collision search O(N log(N))
        } collision;
    /**
     * @brief Available integrators
     */
    enum {
        REB_INTEGRATOR_IAS15 = 0,   ///< IAS15 integrator, 15th order, non-symplectic (default)
        REB_INTEGRATOR_WHFAST = 1,  ///< WHFast integrator, symplectic, 2nd order, up to 11th order correctors
        REB_INTEGRATOR_SEI = 2,     ///< SEI integrator for shearing sheet simulations, symplectic, needs OMEGA variable
        REB_INTEGRATOR_WH = 3,      ///< WH integrator (based on swifter), WHFast is recommended, this integrator is in REBOUND for comparison tests only
        REB_INTEGRATOR_LEAPFROG = 4,    ///< LEAPFROG integrator, simple, 2nd order, symplectic
        REB_INTEGRATOR_HYBRID = 5,  ///< HYBRID Integrator for close encounters (experimental)
        REB_INTEGRATOR_NONE = 6,    ///< Do not integrate anything
        } integrator;

    /**
     * @brief Available boundary conditions
     */
    enum {
        REB_BOUNDARY_NONE = 0,      ///< Do not check for anything (default)
        REB_BOUNDARY_OPEN = 1,      ///< Open boundary conditions. Removes particles if they leave the box 
        REB_BOUNDARY_PERIODIC = 2,  ///< Periodic boundary conditions
        REB_BOUNDARY_SHEAR = 3,     ///< Shear periodic boundary conditions, needs OMEGA variable
        } boundary;

    /**
     * @brief Available gravity routines
     */
    enum {
        REB_GRAVITY_NONE = 0,       ///< Do not calculate graviational forces
        REB_GRAVITY_BASIC = 1,      ///< Basic O(N^2) direct summation algorithm, choose this for shearing sheet and periodic boundary conditions
        REB_GRAVITY_COMPENSATED = 2,    ///< Direct summation algorithm O(N^2) but with compensated summation, slightly slower than BASIC but more accurate
        REB_GRAVITY_TREE = 3,       ///< Use the tree to calculate gravity, O(N log(N)), set opening_angle2 to adjust accuracy.
        } gravity;
    /** @} */


    /**
     * \name Integrator structs (the contain integrator specific variables and temporary data structures) 
     * @{
     */
    struct reb_simulation_integrator_sei ri_sei;        ///< The SEI struct 
    struct reb_simulation_integrator_wh ri_wh;      ///< The WH struct 
    struct reb_simulation_integrator_hybrid ri_hybrid;  ///< The Hybrid struct 
    struct reb_simulation_integrator_whfast ri_whfast;  ///< The WHFast struct 
    struct reb_simulation_integrator_ias15 ri_ias15;    ///< The IAS15 struct 
    /** @} */

    /**
     * \name Callback functions
     * @{
     */
    /**
     * @brief This function allows the user to add additional (non-gravitational) forces.
     */
    void (*additional_forces) (struct reb_simulation* const r);
    /**
     * @brief This function allows the user to modify the dditional (non-gravitational) forces.
     */
    void (*post_timestep_modifications) (struct reb_simulation* const r);
    /**
     * @brief This function is called at the beginning of the simulation and at the end of
     * each timestep.
     */
    void (*heartbeat) (struct reb_simulation* r);
    /**
     * @brief Return the coefficient of restitution. By default it is NULL, assuming a coefficient of 1.
     * @details The velocity of the collision is given to allow for velocity dependent coefficients
     * of restitution.
     */
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v); 

    /**
     * @brief Resolve collision within this function. By default it is NULL, assuming hard sphere model.
     * @details A return value of 0 indicates that both particles remain in the simulation. A return value of 1 (2) indicates that particle 1 (2) should be removed from the simulation. A return value of 3 indicates that both particles should be removed from the simulation. 
     */
    int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);
    /** @} */
    
    /**
     * * \name Hooks for external libraries
     * * @{
     */
    /**
     * @brief Pointer to connect additional (optional) libraries, e.g., reboundx
     */
    void* extras;
    /** @} */
};

/**
 * @name Main REBOUND routines
 * @{
 */
/**
 * @defgroup MainRebFunctions List of the main REBOUND API functions
 * @details These are the functions that typically need to be called by the user.
 * @{
 */
/**
 * @brief Creates and initialises a REBOUND simulation
 * @details Allocate memory for one reb_simulation structure, initialise all variables 
 * and returni the pointer to the reb_simulation sructure. This function must be called 
 * before any particles are added.
 */
struct reb_simulation* reb_create_simulation(void);

/**
 * @brief Initialize reb_simulation structure.
 *
 * @details Same as reb_create_simulation() but does not allocate memory for structure itself.
 * @param r Structure to be initialized (needs to be allocated externally).
 */
void reb_init_simulation(struct reb_simulation* r);

/**
 * @brief Performon one integration step
 * @details You rarely want to call this function yourself.
 * Use reb_integrate instead.
 * @param r The rebound simulation to be integrated by one step.
 */
void reb_step(struct reb_simulation* const r);

/**
 * @brief Performs the actual integration
 * @details This function performs an integration  from the current time t until time tmax.
 * @param r The rebound simulation to be integrated.
 * @param tmax The time to be integrated to. Set this to INFINITY to integrate forever.
 * @return This function returns an integer, indicating the success of the integration.
 */
enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax);

/**
 * @brief Synchronize particles manually at end of timestep
 * @details This function should be called if the WHFAST integrator
 * is used, safe_mode is set to zero and an output is needed.
 * This advances the positions and velocities to be synchronized.
 * If enabled, it also applies the symplectic corrector.
 * If safe_mode is enabled, this function has no effect.
 * @param r The rebound simulation to be synchronized
 */
void reb_integrator_synchronize(struct reb_simulation* r);

/** 
 * @brief Cleanup all temporarily stored integrator values.
 * @param r The rebound simulation to be considered
 **/
void reb_integrator_reset(struct reb_simulation* r);

/**
 * @brief Configure the boundary/root box
 * @details This function helps to setup the variables for the simulation box.
 * Call this function when using open, periodic or shear periodic boundary conditions.
 * @param r The rebound simulation to be considered
 * @param boxsize The size of the root box
 * @param root_nx The numbe rof root boxes in the x direction.
 * @param root_ny The numbe rof root boxes in the y direction.
 * @param root_nz The numbe rof root boxes in the z direction.
 */
void reb_configure_box(struct reb_simulation* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

/**
 * @brief Frees up all space used by a REBOUND simulation and the reb_simulation structure itself.
 * @details The REBOUND simulation is not usable anymore after being passed to this function.
 * @param r The rebound simulation to be freed
 */
void reb_free_simulation(struct reb_simulation* const r);

/**
 * @brief Frees up all space used by a REBOUND simulation, but not the reb_simulation structure itself.
 * @details The REBOUND simulation is not usable anymore after being passed to this function.
 * @param r The rebound simulation to be freed
 */
void reb_free_pointers(struct reb_simulation* const r);

#ifdef MPI
/**
 * @brief Init MPI for simulation r
 */
void reb_mpi_init(struct reb_simulation* const r);

/**
 * @brief Finalize MPI for simulation r
 */
void reb_mpi_finalize(struct reb_simulation* const r);
#endif // MPI

/**
 * @cond PRIVATE
 */
/**
 * @brief Function used to allow binary input.
 */
void reb_reset_temporary_pointers(struct reb_simulation* const r);
/**
 * @brief Function used to allow binary input.
 */
void reb_reset_function_pointers(struct reb_simulation* const r);
/** @endcond */

/** 
 * @brief Adds a particle to the simulation. 
 * @details This function adds the particle pt to the simulation. 
 * @param r The rebound simulation to which the particle will be added
 * @param pt The particle to be added. Note that this is a structure, not a reference to a structure.
 */
void reb_add(struct reb_simulation* const r, struct reb_particle pt);


/**
 * @brief Remove all particles
 * @param r The rebound simulation to be considered
 */
void reb_remove_all(struct reb_simulation* const r);

/**
 * @brief Remove a particle by the position in particles array
 * @param r The rebound simulation to be considered
 * @param index The index in the particles array of the particle to be removed.
 * @param keepSorted Set to 1, then particles with indices higher than index
 * are all shifted down one position, ensuring the ordering remains.
 * @return Returns 1 if particle was successfully removed, 0 if index passed was 
 * out of range.
 */
int reb_remove(struct reb_simulation* const r, int index, int keepSorted);

/**
 * @brief Remove a particle by its id.
 * @param r The rebound simulation to be considered
 * @param id The id of the particle to be removed.
 * @param keepSorted If set to 1 keep the particles with indices in the particles array
 * higher than the one with the passed id are all shifted down one position,
 * ensuring the ordering remains. 
 * @return Returns 1 if particle successfully removed,
 * 0 if id was not found in the particles array.
 */
int reb_remove_by_id(struct reb_simulation* const r, int id, int keepSorted);

/**
 * @brief Run the heartbeat function and check for escaping/colliding particles.
 * @details You rarely want to call this function yourself. It is used internally to 
 * call the function you set to the heartbeat variable in reb_simulation.
 * @param r The rebound simulation to be considered
 */
void reb_run_heartbeat(struct reb_simulation* const r);

/**
 * @brief Hardsphere collision resolving routine (default).
 */
int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);

/**
 * @brief Merging collision resolving routine.
 * @details Merges particle with higher index into particle of lower index.
 *          Conserves mass, momentum and volume. Compatible with HYBARID. 
 */
int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);

/** @} */
/** @} */





/**
 * \name Tools
 * @{
 */
/**
 * @defgroup ToolsRebFunctions List of the helper functions for REBOUND
 * @{
 */
/**
 * @brief Return uniformly distributed random variable in a given range.
 * @param min Minimum value.
 * @param max Maximum value.
 * @return A random variable
 */
double reb_random_uniform(double min, double max);

/**
 * @brief Returns a random variable drawn form a powerlaw distribution.
 * @param min Minimum value.
 * @param max Maximum value.
 * @param slope Slope of powerlaw distribution.
 * @return A random variable
 */
double reb_random_powerlaw(double min, double max, double slope);

/**
 * @brief Return a random number with normal distribution.
 * @details Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
 * @param variance Variance of normal distribution.
 * @return A random variable
 */
double reb_random_normal(double variance);

/**
 * @brief Return a random variable drawn form a Rayleigh distribution.  
 * @details Calculated as described on Rayleigh distribution wikipedia page
 * @param sigma Scale parameter.
 * @return A random variable
 */
double reb_random_rayleigh(double sigma);

/**
 * @brief Move to center of momentum and center of mass frame.
 * @details This function moved all particles to the center of mass 
 * frame (sometimes also called center of momentum frame). In this frame
 * the center of mass is at rest.
 * It is recommended to call this function before you are doing a long
 * term orbit integration. If the particles are slowly drifting away from the
 * coordinate origin, numerical errors might build up.
 * @param r The rebound simulation to be considered
 */
void reb_move_to_com(struct reb_simulation* const r);

/**
 * @brief Returns the center of mass.
 * @param r The rebound simulation to be considered
 * @return The center of mass as a particle (mass, position and velocity correspond to the center of mass)
 */
struct reb_particle reb_get_com(struct reb_simulation* r);

/**
 * @brief Returns the center of mass of two particles
 * @param p1 One of the two particles
 * @param p2 One of the two particles
 * @return The center of mass as a particle (mass, position and velocity correspond to the center of mass)
 */
struct reb_particle reb_get_com_of_pair(struct reb_particle p1, struct reb_particle p2);
/** @} */
/** @} */

/**
 * @brief Takes the center of mass of a system of particles and returns the center of mass with one of the particles removed. 
 * @param com A particle structure that holds the center of mass state for a system of particles (mass, position, velocity).
 * @param p The particle to be removed from com.
 * @return The center of mass with particle p removed.
 */

struct reb_particle reb_get_com_without_particle(struct reb_particle com, struct reb_particle p);

/**
 * @brief Returns a particle pointer's index in the simulation it's in.
 * @param p A pointer to the particle 
 * @return The integer index of the particle in its simulation (will return -1 if not found in the simulation).
 */
int reb_get_particle_index(struct reb_particle* p);

/**
 * @brief Returns the center of mass for particles with indices between first (inclusive) and last (exclusive).
 * @details For example, reb_get_com_range(r, 6, 9) returns COM for particles 6, 7 and 8. 
 * @param r A pointer to the simulation structure.
 * @param first First index in range to consider.
 * @param last Will consider particles with indices < last (i.e., particle with index last not considered).
 * @return A reb_particle structure for the center of mass of all particles in range [first, last). Returns particle filled with zeros if passed last <= first.
 */

struct reb_particle reb_get_com_range(struct reb_simulation* r, int first, int last);

/**
 * @brief Returns the jacobi center of mass for a given particle
 * @param p A pointer to the particle
 * @return A reb_particle structure for the center of mass of all particles with lower index.  Returns particles[0] if passed the 0th particle.
 */

struct reb_particle reb_get_jacobi_com(struct reb_particle* p);

/**
 * \name Built-in output function
 * @{
 */
/**
 * @defgroup OutputRebFunctions List of the built-in output functions for REBOUND
 * @{
 */
/**
 * @brief This function checks if a new output is required at this time.
 * @details This is typically used within the heartbeat function to generate
 * equally spaced outputs.
 * @param interval Output interval.
 * @param r The rebound simulation to be considered
 * @return The return value is 1 if an output is required and 0 otherwise.
 */
int reb_output_check(struct reb_simulation* r, double interval);

/**
 * @brief Output status information on the screen.
 * @details Outputs the current number of particles, the time and the time difference since the last output to the screen.
 * @param r The rebound simulation to be considered
 * @param tmax The maximum integration time (used to calculate the progress in percent)
 */
void reb_output_timing(struct reb_simulation* r, const double tmax);

/**
 * @brief Append an ASCII file with orbital paramters of all particles.
 * @details The orbital parameters are calculated with respect the center of mass.
 * reb_particles are assumed to be sorted from the inside out, the central object having index 0. 
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_orbits(struct reb_simulation* r, char* filename);

/**
 * @brief Save the reb_simualtion structure as a binary
 * @details This function can be used to save the current status of a REBOUND simualtion 
 * and later restart the simualtion.
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_binary(struct reb_simulation* r, char* filename);

/**
 * @brief Append the positions and velocities of all particles to an ASCII file.
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_ascii(struct reb_simulation* r, char* filename);

/**
 * @brief Write the positions of all particles to a binary file.
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_binary_positions(struct reb_simulation* r, char* filename);

/**
 * @brief Append the velocity dispersion of the particles to an ASCII file.
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename);
/** @} */
/** @} */



/**
 * \name Built-in setup/input functions
 * @{
 */
/**
 * @defgroup SetupRebFunctions List of the built-in setup helper functions for REBOUND
 * @{
 */
/**
 * @brief returns the true anomaly for a given eccentricity and mean anomaly
 * @param e Eccentricity
 * @param M Mean anomaly
 * @return True anomaly
 */
double reb_tools_M_to_f(double e, double M);

/**
 * @brief Initialize a particle on an orbit in the xy plane.
 * @param G Gravitational constant.
 * @param primary Particle structure for the orbit's reference body.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param omega Pericenter of the particle.
 * @param f true anomaly of the particle.
 * @return Returns a particle structure with the given orbital parameters. 
 */
struct reb_particle reb_tools_orbit2d_to_particle(double G, struct reb_particle primary, double m, double a, double e, double omega, double f);

/**
 * @brief Initialize a particle on a 3D orbit, passing an error variable to flag why particle is set to nan.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @details Error codes:\n
 * 1. Can't set e exactly to 1.\n
 * 2. Eccentricity can never be less than zero.\n
 * 3. Bound orbit (a>0) can't have e>1.\n
 * 4. Unbound orbit (a<0) can't have e<1.\n
 * 5. Unbound orbit can't have f set beyond the asymptotes defining the particle.
 * @param G Gravitational constant.
 * @param primary Particle structure for the orbit's reference body.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane..
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 * @param err Pointer to error code that wil be set by this function. Used for checking why particle was set to nans.
 * @return Returns a particle structure with the given orbital parameters. 
 */
struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f, int* err);

/**
 * @brief Initialize a particle on a 3D orbit.  See Fig. 2.13 of Murray & Dermott Solar System Dynamics for diagram.
 * @param G Gravitational constant.
 * @param primary Particle structure for the orbit's reference body.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param e Eccentricity of the particle.
 * @param i inclination of the particle to the reference plane.
 * @param Omega Longitude of the ascending node of the particle.
 * @param omega argument of pericenter of the particle.
 * @param f true anomaly of the particle.
 * @return Returns a particle structure with the given orbital parameters. 
 */
struct reb_particle reb_tools_orbit_to_particle(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f);

/**
 * @brief This function calculates orbital elements for a given particle, passing an error variable to flag why orbit is set to nan.
 * @details Error codes:\n
 * 1. Primary has no mass.\n
 * 2. Particle and primary positions are the same.\n
 * @param G The gravitational constant.
 * @param p reb_particle for which the orbit is calculated.
 * @param primary Particle structure for the orbit's reference body.
 * @param err error code for checking why orbit was set to nans.
 * @return reb_orbit struct with orbital parameters. 
 */
struct reb_orbit reb_tools_particle_to_orbit_err(double G, struct reb_particle p, struct reb_particle primary, int* err);

/**
 * @brief This function calculates orbital elements for a given particle. 
 * @param G The gravitational constant.
 * @param p reb_particle for which the orbit is calculated.
 * @param primary Particle structure for the orbit's reference body.
 * @return reb_orbit struct with orbital parameters. 
 */
struct reb_orbit reb_tools_particle_to_orbit(double G, struct reb_particle p, struct reb_particle primary);

/**
 * @brief Initialize a particle on a 3D orbit.  See Pal 2009 for a definition of these coordinates.
 * @detail Pal describes a coordinate system for Keplerian Orbits that is analytical (i.e. infinitely differentiable) between spatial coordinates and orbital elements. See http://adsabs.harvard.edu/abs/2009MNRAS.396.1737P
 * @param G Gravitational constant.
 * @param primary Particle structure for the orbit's reference body.
 * @param m Mass of the particle.
 * @param a Semi-major axis of the particle.
 * @param lambda longitude.
 * @param k Eccentricity/pericenter k = e*cos(w).
 * @param h Eccentricity/pericenter h = e*sin(w).
 * @param ix Inclination, x component.
 * @param iy Inclination, y component.
 * @return Returns a particle structure with the given orbital parameters. 
 */
struct reb_particle reb_tools_pal_to_particle(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy);


/**
 * @brief Reads a binary file.
 * @details Also initialises the particles array with data form the binary file.
 * This can be used to restart a simualtion.
 * @param filename Filename to be read.
 * @return Returns a pointer to a REBOUND simulation.
 */
struct reb_simulation* reb_create_simulation_from_binary(char* filename);

/**
 * @brief This function sets up a Plummer sphere.
 * @param r The rebound simulation to be considered
 * @param _N Number of particles in the plummer sphere.
 * @param M Total mass of the cluster.
 * @param R Characteristic radius of the cluster.
 */
void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R);

/**
 * @brief Reads arguments from the command line.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @param argument Argument to look for.
 * @return Returns NULL if argument was not given. Return the argument otherwise.
 */
char* reb_read_char(int argc, char** argv, const char* argument);

/**
 * @brief Reads arguments as a double value from the command line.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @param argument Argument to look for.
 * @param _default Default value.
 * @return Returns _default if argument was not given. Return the argument converted to double otherwise.
 */
double reb_read_double(int argc, char** argv, const char* argument, double _default);

/**
 * @brief Reads arguments as a int value from the command line.
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @param argument Argument to look for.
 * @param _default Default value.
 * @return Returns _default if argument was not given. Return the argument converted to int otherwise.
 */
int reb_read_int(int argc, char** argv, const char* argument, int _default);
/** @} */
/** @} */



/**
 * \name Setup functions for variational particles.
 * @{
 * @brief This function calculates the first/second derivative of a Keplerian orbit. 
 * @details Derivatives of Keplerian orbits are required for variational equations, in particular
 *          for optimization problems. 
 *          The derivative is calculated with respect to the variables that appear in the function name.
 *          One variable implies that a first derivative is returned, two variables implies that a second
 *          derivate is returned. Classical orbital parameters and those introduced by Pal (2009) are 
 *          supported. Pal coordinates have the advantage of being analytical (i.e. infinite differentiable).
 *          Classical orbital parameters may have singularities, for example when e is close to 0.
 *          Note that derivatives with respect to Cartesian coordinates are trivial and therefore not
 *          implemented as seperate functions. 
 *          The following variables are supported: a, e, inc, f, omega, Omega, h, k, ix, iy and m (mass). 
 * @return The derivative as a particle structre. Each structure element is a derivative.
 * @param G The gravitational constant
 * @param primary The primary of the Keplerian orbit
 * @param po The original partical for which the derivative is to be calculated.
 */
struct reb_particle reb_derivatives_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_h(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k_k(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_h_h(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_lambda_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_h_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k_h(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_a(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_ix_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_iy_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_h_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_lambda_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_lambda_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_h_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_k_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_ix_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_h(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_k(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_a(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_lambda(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_h(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_k(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_ix(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_iy(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_m(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e_e(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_inc(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_inc_inc(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_Omega_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_omega_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_f_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_e(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_inc(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_a_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e_inc(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_e_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_e(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_inc_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_inc_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_inc_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_inc(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_omega_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_Omega_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_Omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_omega_f(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_omega(double G, struct reb_particle primary, struct reb_particle po);
struct reb_particle reb_derivatives_m_f(double G, struct reb_particle primary, struct reb_particle po);
/** @} */

/**
 * \name Particle manipulation functions
 * @{
 */
/**
 * @defgroup ParticleManipFunctions List of reb_particle manipulation functions for REBOUND
 * @{
 */
/**
 * @brief Subtract particle p2 from particle p1 (p1 - p2).
 * @details Subtracts positions, velocities, accelerations and mass element by element. 
 * @param p1 First reb_particle.
 * @param p2 Second reb_particle to subtract from p1.
 * @returns A new particle with no pointers (not in any simulation etc.) set.
 */
struct reb_particle reb_particle_minus(struct reb_particle p1, struct reb_particle p2);

/**
 * @brief Add particle p1 to particle p1.
 * @details Adds positions, velocities, accelerations and mass element by element. 
 * @param p1 First reb_particle.
 * @param p2 Second reb_particle.
 * @returns A new particle with no pointers (not in any simulation etc.) set.
 */
struct reb_particle reb_particle_plus(struct reb_particle p1, struct reb_particle p2);

/**
 * @brief Multiply a particle's members by a constant.
 * @brief Multiplies particle's positions, velocities, accelerations and mass by a constant.
 * @param p1 reb_particle to modify.
 * @param value Value by which to multiply particle's fields.
 * @returns A new particle with no pointers (not in any simulation etc.) set.
 */
struct reb_particle reb_particle_multiply(struct reb_particle p1, double value);

/**
 * @brief Divide a particle's members by a constant.
 * @brief Divides particle's positions, velocities, accelerations and mass by a constant.
 * @param p1 reb_particle to modify.
 * @param value Value by which to divide particle's fields.
 * @returns A new particle with no pointers (not in any simulation etc.) set.
 */
struct reb_particle reb_particle_divide(struct reb_particle p1, double value);
/** @} */
/** @} */

/**
 * \name Miscellaneous tools
 * @{
 */
/**
 * @defgroup MiscRebFunctions List of the miscellaneous helper functions for REBOUND
 * @{
 */
/**
 * @brief Calculate the total energy (potential and kinetic).
 * @details Does not work for WH/SEI.
 * @param r The rebound simulation to be considered.
 * @return Total energy. 
 */
double reb_tools_energy(const struct reb_simulation* const r);

/**
 * @brief Calculate the system's angular momentum.
 * @param r The rebound simulation to be considered.
 * @return The angular momentum vector as a reb_vec3d struct.
 */
struct reb_vec3d reb_tools_angular_momentum(const struct reb_simulation* const r);

/**
 * @brief Add and initialize a set of first order variational particles
 * @param r The rebound simulation to be considered
 * @param testparticle This flag determines if the set of variational particles is for a testparticle or not.
 * If testparticle is >= 0, then only one variational particle (the test particle) will be added.
 * If testparticle is -1, one variational particle for each real particle will be added.
 * @return Returns the index of the first variational particle added
 **/
int reb_add_var_1st_order(struct reb_simulation* const r, int testparticle);

/**
 * @brief Add and initialize a set of second order variational particles
 * @details Note that a set of second order variational particles requires two sets of first order variational equations.
 * @param r The rebound simulation to be considered
 * @param testparticle This flag determines if the set of variational particles is for a testparticle or not.
 * If testparticle is >= 0, then only one variational particle (the test particle) will be added.
 * If testparticle is -1, one variational particle for each real particle will be added.
 * @param index_1st_order_a The index of the corresponding first variational particles.
 * @param index_1st_order_b The index of the corresponding first variational particles.
 * @return Returns the index of the first variational particle added
 **/
int reb_add_var_2nd_order(struct reb_simulation* const r, int testparticle, int index_1st_order_a, int index_1st_order_b);

/** 
 * @brief Init the MEGNO particles, enable MEGNO calculation
 * @param r The rebound simulation to be considered
 */
void reb_tools_megno_init(struct reb_simulation* const r);

/**
 * @brief Get the current MEGNO value
 * @param r The rebound simulation to be considered
 * @return Returns the current value of the MEGNO
 */
double reb_tools_calculate_megno(struct reb_simulation* r);

/**
 * @brief Returns the largest Lyapunov characteristic number (LCN), or maximal Lyapunov exponent
 * @details MEGNO needs to be enabled to calculate this value.
 * @param r The rebound simulation to be considered
 * @return Returns the current CN
 */
double reb_tools_calculate_lyapunov(struct reb_simulation* r);

/**
 * @brief Print out an error message, then exit in a semi-nice way.
 */
void reb_exit(const char* const msg);

/**
 * @brief Print out a warningr message, then continue.
 */
void reb_warning(const char* const msg);
/** @} */
/** @} */

#endif
