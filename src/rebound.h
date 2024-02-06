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

#ifdef _WIN64
#define _LP64
#endif
#ifdef _WIN32
#include <WinSock2.h>
#define _WINSOCKAPI_ //stops windows.h including winsock.h
#include <windows.h>
#define REB_RESTRICT
#define DLLEXPORT __declspec(dllexport)
#define __restrict__
#elif __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#define REB_RESTRICT
#define DLLEXPORT
#else // Linux and MacOS
#define REB_RESTRICT restrict
#define DLLEXPORT
#endif // _WIN32

#include <stdio.h>
#include <inttypes.h>
#include <stdint.h>
#include <signal.h>
#define _USE_MATH_DEFINES
#include <math.h>

#ifdef _WIN32
typedef struct reb_timeval {
    int64_t tv_sec;
    int64_t tv_usec;
} reb_timeval;
int gettimeofday(struct reb_timeval * tp, struct timezone * tzp);
int asprintf(char **strp, const char *fmt, ...);
int rand_r (unsigned int *seed);
#include <io.h>
#define _TIMEVAL_DEFINED
#else // Linux and MacOS
#define reb_timeval timeval
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#endif // _WIN32

#ifdef AVX512
#include <immintrin.h>
#endif

#ifdef MPI
#include "mpi.h"
#endif

#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif


#ifndef __GNUC__
#  define  __attribute__(x)  /*Ignore attributes in non-GNU compilers*/
#endif


// Global constants and variables
DLLEXPORT extern const char* reb_build_str;   ///< Date and time build string.
DLLEXPORT extern const char* reb_version_str; ///< Version string.
DLLEXPORT extern const char* reb_githash_str; ///< Current git hash.
DLLEXPORT extern const char* reb_logo[26];    ///< Logo of rebound. 
DLLEXPORT extern const unsigned char reb_favicon_png[]; /// < Favicon in PNG format.
DLLEXPORT extern const unsigned int reb_favicon_len;
DLLEXPORT extern const int reb_max_messages_length;
DLLEXPORT extern const int reb_N_max_messages;
extern volatile sig_atomic_t reb_sigint;  ///< Graceful global interrupt handler 

// Forward declarations
struct reb_simulation;
struct reb_simulationarchive;
struct reb_display_data;
struct reb_server_data;
struct reb_treecell;
struct reb_variational_configuration;
struct reb_display_settings;

// Particle structure
struct reb_particle {
    double x;                   // Cartesian coordinates
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    double m;                   // mass
    double r;                   // physical radius
    double last_collision;      // Last time the particle had a physical collision.
    struct reb_treecell* c;     // Pointer to the cell the particle is currently in.
#if !defined(_LP64)
    char pad1[4];   // c is short by 4 bytes
#endif
    uint32_t hash;              // Hash, can be used to identify particle.
#if !defined(_LP64)
    char pad2[4];   // ap is not padded to 8 bytes
#endif
    void* ap;                   // This pointer allows REBOUNDx to add additional properties to the particle.
#if !defined(_LP64)
    char pad3[4];   // ap is short by 4 bytes
#endif
    struct reb_simulation* sim; // Pointer to the parent simulation.
#if !defined(_LP64)
    char pad4[4];   // sim is short by 4 bytes
#endif
};

// Generic 3d vector
struct reb_vec3d {
    double x;
    double y;
    double z;
};

// Generic 4d matrix (single precision)
struct reb_mat4df {
    float m[16];
};


// Generic 6d vector
struct reb_vec6d{
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

// Rotation (implemented as a quaternion)
struct reb_rotation {
    double ix;
    double iy;
    double iz;
    double r;
};

// Structure representing one particle-particle collision
struct reb_collision{
    int p1;                 // Index of first particle involved in collision
    int p2;                 // Index of second particle
    struct reb_vec6d gb;    // Offset due to boundary conditions
    int ri;                 // Root cell index (MPI only)
};

// Generic pointer with 7 elements, for internal use only (IAS15).
struct reb_dp7 {
    double* REB_RESTRICT p0;
    double* REB_RESTRICT p1;
    double* REB_RESTRICT p2;
    double* REB_RESTRICT p3;
    double* REB_RESTRICT p4;
    double* REB_RESTRICT p5;
    double* REB_RESTRICT p6;
};

// Integrator structures 
// IAS15 (Rein & Spiegel 2015)
struct reb_integrator_ias15 {
    double epsilon;                         // Precision control parameter
    double min_dt;                          // Minimal timestep
    unsigned int adaptive_mode;             // 0: fractional error is calculated seperately for each particle
                                            // 1: fractional error is calculated globally (default)
                                            // 2: Dang, Rein & Spiegel (2023) timestep criterion
                                            // 3: Aarseth (1985) timestep criterion
    // Internal use
    uint64_t iterations_max_exceeded; // Counter how many times the iteration did not converge. 
    unsigned int N_allocated;          
    double* REB_RESTRICT at;
    double* REB_RESTRICT x0;
    double* REB_RESTRICT v0;
    double* REB_RESTRICT a0;
    double* REB_RESTRICT csx;
    double* REB_RESTRICT csv;
    double* REB_RESTRICT csa0;
    struct reb_dp7 g;
    struct reb_dp7 b;
    struct reb_dp7 csb;             // Compensated summation storage for b
    struct reb_dp7 e;
    struct reb_dp7 br;              // Used for resetting the b coefficients if a timestep gets rejected
    struct reb_dp7 er;              // Same for e coefficients
    int* map;                       // internal map to particles (this is an identity map except when MERCURIUS is used
    unsigned int N_allocated_map;   // allocated size for map
};

// Mercurius (Rein et al. 2019)
struct reb_integrator_mercurius {
    double (*L) (const struct reb_simulation* const r, double d, double dcrit); // Switching function (default same as Mercury) 
    double r_crit_hill;                                 // Critical switching distance in units of Hill radii
    unsigned int recalculate_coordinates_this_timestep; // Set to 1 if particles have been modified
    unsigned int recalculate_r_crit_this_timestep;      // Set to 1 if to recalculate critical switching radii
    unsigned int safe_mode;                             // Combine Kick steps at beginning and end of timestep
   
    // Internal use
    unsigned int is_synchronized;   
    unsigned int mode;              // 0 if WH is operating, 1 if IAS15 is operating.
    unsigned int encounter_N;       // Number of particles currently having an encounter
    unsigned int encounter_N_active;// Number of active particles currently having an encounter
    unsigned int tponly_encounter;  // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles
    unsigned int N_allocated;
    unsigned int N_allocated_additional_forces;
    unsigned int N_allocated_dcrit; // Current size of dcrit arrays
    double* dcrit;                  // Precalculated switching radii for particles
    struct reb_particle* REB_RESTRICT particles_backup; //  contains coordinates before Kepler step for encounter prediction
    struct reb_particle* REB_RESTRICT particles_backup_additional_forces; // contains coordinates before Kepler step for encounter prediction
    int* encounter_map;             // Map to represent which particles are integrated with ias15
    struct reb_vec3d com_pos;       // Used to keep track of the center of mass during the timestep
    struct reb_vec3d com_vel;
};

// Symplectic Epicycle Integrator SEI (Rein & Tremaine 2011)
struct reb_integrator_sei {
    double OMEGA;       // Epicyclic frequency
    double OMEGAZ;      // Epicyclic frequency in z direction (if not set, use OMEGA)

    // Internal use
    double lastdt;      // Cached sin(), tan() for this value of dt.
    double sindt;       // Cached sin() 
    double tandt;       // Cached tan() 
    double sindtz;      // Cached sin(), z axis
    double tandtz;      // Cached tan(), z axis
};

// SABA Integrator (Laskar & Robutel 2001)
struct reb_integrator_saba {
    enum {
        REB_SABA_1 = 0x0, // WH
        REB_SABA_2 = 0x1, // SABA2
        REB_SABA_3 = 0x2, // SABA3
        REB_SABA_4 = 0x3, // SABA4
        REB_SABA_CM_1 = 0x100, // SABACM1 (Modified kick corrector)
        REB_SABA_CM_2 = 0x101, // SABACM2 (Modified kick corrector)
        REB_SABA_CM_3 = 0x102, // SABACM3 (Modified kick corrector)
        REB_SABA_CM_4 = 0x103, // SABACM4 (Modified kick corrector)
        REB_SABA_CL_1 = 0x200, // SABACL1 (lazy corrector)
        REB_SABA_CL_2 = 0x201, // SABACL2 (lazy corrector)
        REB_SABA_CL_3 = 0x202, // SABACL3 (lazy corrector)
        REB_SABA_CL_4 = 0x203, // SABACL4 (lazy corrector)
        REB_SABA_10_4 = 0x4,   // SABA(10,4), 7 stages
        REB_SABA_8_6_4 = 0x5,  // SABA(8,6,4), 7 stages
        REB_SABA_10_6_4 = 0x6, // SABA(10,6,4), 8 stages, default
        REB_SABA_H_8_4_4 = 0x7,// SABAH(8,4,4), 6 stages
        REB_SABA_H_8_6_4 = 0x8,// SABAH(8,6,4), 8 stages
        REB_SABA_H_10_6_4 = 0x9,// SABAH(10,6,4), 9 stages
    } type;                             // Type of integrator
    unsigned int safe_mode;             // Combine first and last sub-step
    unsigned int is_synchronized;       // 1: physical state, 0: needs synchronization
    unsigned int keep_unsynchronized;   // 1: continue from unsynchronized state after synchronization
};

// WHFast Integrator (Rein & Tamayo 2015)
struct reb_integrator_whfast {
    unsigned int corrector;                                     // Order of first symplectic corrector: 0 (default - no corrector), 3, 5, 7, 11, 17.  
    unsigned int corrector2;                                    // 0: no second corrector, 1: use second corrector
    enum {
        REB_WHFAST_KERNEL_DEFAULT = 0,
        REB_WHFAST_KERNEL_MODIFIEDKICK = 1,
        REB_WHFAST_KERNEL_COMPOSITION = 2,
        REB_WHFAST_KERNEL_LAZY = 3,
    } kernel;                                                   // Kernel type. See Rein, Tamayo & Brown 2019 for details.
    enum {
        REB_WHFAST_COORDINATES_JACOBI = 0,                      // Jacobi coordinates (default)
        REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC = 1,      // Democratic Heliocentric coordinates
        REB_WHFAST_COORDINATES_WHDS = 2,                        // WHDS coordinates (Hernandez and Dehnen, 2017)
    } coordinates;                                              // Coordinate system used in Hamiltonian splitting
    unsigned int recalculate_coordinates_this_timestep;         // 1: recalculate coordinates from inertial coordinates
    unsigned int safe_mode;                                     // 0: Drift Kick Drift scheme (default), 1: combine first and last sub-step.
    unsigned int keep_unsynchronized;                           // 1: continue from unsynchronized state after synchronization 

    // Internal use
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    struct reb_particle* REB_RESTRICT p_temp;   // Used for lazy implementer's kernel 
    unsigned int is_synchronized;
    unsigned int N_allocated;
    unsigned int N_allocated_tmp;               // Used for lazy implementer's kernel 
    unsigned int timestep_warning;
    unsigned int recalculate_coordinates_but_not_synchronized_warning;
};

// Special particle struct for WHFast512
struct reb_particle_avx512{
#ifdef AVX512
    __m512d m __attribute__ ((aligned (64)));
    __m512d x __attribute__ ((aligned (64)));
    __m512d y __attribute__ ((aligned (64)));
    __m512d z __attribute__ ((aligned (64)));
    __m512d vx __attribute__ ((aligned (64)));
    __m512d vy __attribute__ ((aligned (64)));
    __m512d vz __attribute__ ((aligned (64)));
#else // AVX512
    double m[8]; // dummy for when AVX512 is not available
    double x[8];
    double y[8];
    double z[8];
    double vx[8];
    double vy[8];
    double vz[8];
#endif // AVX512
};

// WHFast512 Integrator (Javaheri & Rein 2023)
struct reb_integrator_whfast512 {
    unsigned int gr_potential;          // 1: Turn on GR potential of central object, 0 (default): no GR potential
    unsigned int N_systems;             // Number of systems to be integrator in parallel: 1 (default, up to 8 planets), 2 (up to 4 planets each), 4 (2 planets each)
    unsigned int keep_unsynchronized;   // 1: continue from unsynchronized state after synchronization 

    // Internal use
    unsigned int is_synchronized;
    unsigned int N_allocated;
    unsigned int recalculate_constants;
    struct reb_particle_avx512* p_jh;
    struct reb_particle p_jh0[4];
};

// Bulirsch Stoer Integrator (roughly follows fortran code by E. Hairer and G. Wanner)
struct reb_integrator_bs {
    double eps_abs; // Allowed absolute scalar error.
    double eps_rel; // Allowed relative scalar error.
    double min_dt;  // Minimum timestep
    double max_dt;  // Maximum teimstep

    // Internal use
    struct reb_ode* nbody_ode;  // ODE corresponding to N-body system
    int* sequence;              // stepsize sequence
    int* cost_per_step;         // overall cost of applying step reduction up to iteration k + 1, in number of calls.
    double* cost_per_time_unit; // cost per unit step.
    double* optimal_step;       // optimal steps for each order. 
    double* coeff;              // extrapolation coefficients.
    double dt_proposed;
    int first_or_last_step;
    int previous_rejected;
    int target_iter;
    int user_ode_needs_nbody;   // Do not set manually. Use needs_nbody in reb_ode instead.
};

// Available methods for EOS Integrator
enum REB_EOS_TYPE {
    REB_EOS_LF = 0x00, 
    REB_EOS_LF4 = 0x01,
    REB_EOS_LF6 = 0x02,
    REB_EOS_LF8 = 0x03, 
    REB_EOS_LF4_2 = 0x04,
    REB_EOS_LF8_6_4= 0x05,
    REB_EOS_PLF7_6_4= 0x06,
    REB_EOS_PMLF4 = 0x07,
    REB_EOS_PMLF6 = 0x08,
};

// Embedded Operator Splitting Integrator (Rein 2020)
struct reb_integrator_eos {
    enum REB_EOS_TYPE phi0;         // Outer operator splitting method
    enum REB_EOS_TYPE phi1;         // Inner operator splitting method
    unsigned int n;                 // Number of inner splittings per outer splitting
    unsigned int safe_mode;         // Combine Kick steps at beginning and end of timestep

    // Internal use
    unsigned int is_synchronized;
};


// Integer-based positions and velocities for particles. Used in JANUS integrator. 
#define REB_PARTICLE_INT_TYPE int64_t
struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};

// Janus integrator (Rein & Tamayo 2018)
struct reb_integrator_janus {
    double scale_pos;       // Scale of position grid. Default 1e-16
    double scale_vel;       // Scale of velocity grid. Default 1e-16
    unsigned int order;     // Order: 2 (default), 4, 6, 8, 10 
    unsigned int recalculate_integer_coordinates_this_timestep;  // Set to 1 if particles have been modified

    // Internal use
    struct reb_particle_int* REB_RESTRICT p_int;
    unsigned int N_allocated;
};

// Possible return values of of rebound_integrate
enum REB_STATUS {
    // Any status less than SINGLE_STEP get incremented once every timestep until SINGLE_STEP is reached.
    REB_STATUS_SINGLE_STEP = -10, // Performing a single step, then switching to PAUSED.
    REB_STATUS_PAUSED = -3,       // Simulation is paused by visualization.
    REB_STATUS_LAST_STEP = -2,    // Current timestep is the last one. Needed to ensure that t=tmax exactly.
    REB_STATUS_RUNNING = -1,      // Simulation is current running, no error occurred.
    REB_STATUS_SUCCESS = 0,       // Integration finished successfully.
    REB_STATUS_GENERIC_ERROR = 1, // A generic error occurred and the integration was not successful.
    REB_STATUS_NO_PARTICLES = 2,  // The integration ends early because no particles are left in the simulation.
    REB_STATUS_ENCOUNTER = 3,     // The integration ends early because two particles had a close encounter (see exit_min_distance)
    REB_STATUS_ESCAPE = 4,        // The integration ends early because a particle escaped (see exit_max_distance)  
    REB_STATUS_USER = 5,          // User caused exit, simulation did not finish successfully.
    REB_STATUS_SIGINT = 6,        // SIGINT received. Simulation stopped.
    REB_STATUS_COLLISION = 7,     // The integration ends early because two particles collided. 
};

// Holds a particle's hash and the particle's index in the particles array. Used for particle_lookup_table.
struct reb_hash_pointer_pair{
    uint32_t hash;
    int index;
};

// Main REBOUND Simulation structure
// Note: only variables that should be accessed by users are documented here.
struct reb_simulation {
    double  t;                      // Current simulation time. Default: 0.
    double  G;                      // Gravitational constant. Default: 1.
    double  softening;              // Gravitational softening. Default: 0.
    double  dt;                     // Timestep. Default: 0.001.
    double  dt_last_done;           // Last successful timestep.
    uint64_t steps_done;  // Number of timesteps done.
    unsigned int     N;             // Number of particles (includes variational particles). Default: 0.
    int     N_var;                  // Number of variational particles. Default 0.
    unsigned int     N_var_config; 
    struct reb_variational_configuration* var_config;   // Configuration structs. These contain details on variational particles. 
    int     var_rescale_warning;    
    int     N_active;               // Number of active (i.e. not test-particle) particles. Default: -1 (all particles are active). 
    int     testparticle_type;      // 0 (default): active particles do not feel test-particles, 1: active particles feel test-particles
    int     testparticle_hidewarnings;
    struct  reb_hash_pointer_pair* particle_lookup_table;
    int     hash_ctr;
    int     N_lookup;               // Number of entries in particle_lookup_table.
    int     N_allocated_lookup;     // Number of lookup table entries allocated.
    unsigned int   N_allocated;     // Current maximum space allocated in the particles array on this node. 
    struct reb_particle* particles; // Main particle array with active, variational, and test particles.
    struct reb_vec3d* gravity_cs; 
    int     N_allocated_gravity_cs;
    struct reb_treecell** tree_root;
    int     tree_needs_update;      // Flag to force a tree update (after boundary check)
    double opening_angle2;          // Opening angle for tree-based gravity calculation. Defaukt 0.25. 
    enum REB_STATUS status;         // Current simulation status
    int     exact_finish_time;      // 1 (default): integrate exactly to the time requested and adjust timestep if needed, 0: may overshoot by one timestep

    unsigned int force_is_velocity_dependent; // 0 (default): force only depends on position, 1: force also depends on velocities
    unsigned int gravity_ignore_terms;
    double output_timing_last;      // Time when reb_simulation_output_timing() was called the last time. 
    int save_messages;              // 0 (default): print messages on screen, 1: ignore messages (used in python interface).
    char** messages;                // Array of strings containing last messages (only used if save_messages==1). 
    double exit_max_distance;       // Exit simulation if a particle is this far away from the origin.
    double exit_min_distance;       // Exit simulation if two particles come this close to each other.
    double usleep;                  // Artificially slow down simulations by this many microseconds each timestep.
    struct reb_display_settings* display_settings;// Optional. Will overwrite settings for visualization. If NULL, UI will determine settings.
    struct reb_display_data* display_data; // Datastructure stores visualization related data. Does not have to be modified by the user. 
    struct reb_server_data* server_data; // Datastructure stores server related data. Does not have to be modified by the user. 
    int track_energy_offset;        // 0 (default): do not track energy offset due to merging/lost particles, 1: track offset
    double energy_offset;           // Only used if track_energy_offset = 1
    double walltime;                // Cumulative walltime of entire integration.
    double walltime_last_step;      // Wall time of last step.
    double walltime_last_steps;     // Average wall time of last step (updated every 0.1s).
    double walltime_last_steps_sum;
    int walltime_last_steps_N;
    uint32_t python_unit_l;         // Only used for when working with units in python.
    uint32_t python_unit_m;         // Only used for when working with units in python.
    uint32_t python_unit_t;         // Only used for when working with units in python.
    
    // Simulation domain and ghost boxes 
    struct  reb_vec3d boxsize;      // Size of the entire simulation box, root_x*boxsize. Set in box_init().
    double  boxsize_max;            // Maximum size of the entire box in any direction. Set in box_init().
    double  root_size;              // Size of a root box. 
    int     N_root;                 // Total number of root boxes in all directions, N_root_x*N_root_y*N_root_z. Default: 1. Set in box_init().
    int     N_root_x;               // Number of ghost boxes in x direction. Do not change manually.
    int     N_root_y;
    int     N_root_z;
    int     N_ghost_x;              // Number of ghost boxes in x direction.
    int     N_ghost_y;
    int     N_ghost_z;

    // MPI Parallelization
#ifdef MPI
    int    mpi_id;                              // Unique id of this node (starting at 0). Used for MPI only.
    int    mpi_num;                             // Number of MPI nodes. Used for MPI only.
    struct reb_particle** particles_send;       // Send buffer for particles. There is one buffer per node. 
    int*   N_particles_send;                    // Current length of particle send buffer. 
    int*   N_particles_send_max;                // Maximal length of particle send beffer before realloc() is needed. 
    struct reb_particle** particles_recv;       // Receive buffer for particles. There is one buffer per node. 
    int*   N_particles_recv;                    // Current length of particle receive buffer. 
    int*   N_particles_recv_max;                // Maximal length of particle receive beffer before realloc() is needed. */

    struct reb_treecell** tree_essential_send;  // Send buffer for cells. There is one buffer per node. 
    int*   N_tree_essential_send;               // Current length of cell send buffer. 
    int*   N_tree_essential_send_max;           // Maximal length of cell send beffer before realloc() is needed. 
    struct reb_treecell** tree_essential_recv;  // Receive buffer for cells. There is one buffer per node. 
    int*   N_tree_essential_recv;               // Current length of cell receive buffer. 
    int*   N_tree_essential_recv_max;           // Maximal length of cell receive beffer before realloc() is needed. 
#endif // MPI

    int collision_resolve_keep_sorted;      // 0 (default): may reorder particles during collisions, 1: keep particles sorted.
    struct reb_collision* collisions;       // Array of current collisions. Do not change manually
    int N_allocated_collisions;
    double minimum_collision_velocity;      // Ensure relative velocity during collisions is at least this much (to avoid particles sinking into each other)
    double collisions_plog;                 // Keeping track of momentum transfer in collisions (for ring simulations)
    double max_radius0;                     // The largest particle radius, set automatically, needed for collision search.
    double max_radius1;                     // The second largest particle radius, set automatically, needed for collision search.
    int64_t collisions_log_n;                  // Cumulative number of collisions in entire simulation.
    
    // MEGNO Chaos indicator. These variables should not be accessed directly. Use functions provided instead.
    int calculate_megno;    // Do not change manually. Internal flag that determines if megno is calculated (default=0, but megno_init() sets it to the index of variational particles used for megno)
    double megno_Ys;        // Running megno sum (internal use)
    double megno_Yss;       // Running megno sum (internal use)
    double megno_cov_Yt;    // covariance of MEGNO Y and t
    double megno_var_t;     // variance of t 
    double megno_mean_t;    // mean of t
    double megno_mean_Y;    // mean of MEGNO Y
    int64_t   megno_n;         // number of covariance updates

    unsigned int rand_seed; // seed for random number generator, used by MEGNO and other random number generators in REBOUND.
    
     // Simulationarchive. These variables should not be accessed directly. Use functions provided instead. 
    int    simulationarchive_version;               // Version of the SA binary format (1=original/, 2=incremental)
    double simulationarchive_auto_interval;         // Current sampling cadence, in code units
    double simulationarchive_auto_walltime;         // Current sampling cadence, in wall time
    uint64_t simulationarchive_auto_step; // Current sampling cadence, in time steps
    double simulationarchive_next;                  // Next output time (simulation tim or wall time, depending on wether auto_interval or auto_walltime is set)
    uint64_t simulationarchive_next_step; // Next output step (only used if auto_steps is set)
    char*  simulationarchive_filename;              // Name of output file

    // Available modules in REBOUND
    enum {
        REB_COLLISION_NONE = 0,         // Do not search for collisions (default)
        REB_COLLISION_DIRECT = 1,       // Direct collision search O(N^2)
        REB_COLLISION_TREE = 2,         // Tree based collision search O(N log(N))
        REB_COLLISION_LINE = 4,         // Direct collision search O(N^2), looks for collisions by assuming a linear path over the last timestep
        REB_COLLISION_LINETREE = 5,     // Tree-based collision search O(N log(N)), looks for collisions by assuming a linear path over the last timestep
        } collision;
    enum {
        REB_INTEGRATOR_IAS15 = 0,       // IAS15 integrator, 15th order, non-symplectic (default)
        REB_INTEGRATOR_WHFAST = 1,      // WHFast integrator, symplectic, 2nd order, up to 11th order correctors
        REB_INTEGRATOR_SEI = 2,         // SEI integrator for shearing sheet simulations, symplectic, needs OMEGA variable
        REB_INTEGRATOR_LEAPFROG = 4,    // LEAPFROG integrator, simple, 2nd order, symplectic
        REB_INTEGRATOR_NONE = 7,        // Do not integrate anything
        REB_INTEGRATOR_JANUS = 8,       // Bit-wise reversible JANUS integrator.
        REB_INTEGRATOR_MERCURIUS = 9,   // MERCURIUS integrator 
        REB_INTEGRATOR_SABA = 10,       // SABA integrator family (Laskar and Robutel 2001)
        REB_INTEGRATOR_EOS = 11,        // Embedded Operator Splitting (EOS) integrator family (Rein 2019)
        REB_INTEGRATOR_BS = 12,         // Gragg-Bulirsch-Stoer 
        // REB_INTEGRATOR_TES = 20,     // Used to be Terrestrial Exoplanet Simulator (TES) -- Do not reuse.
        REB_INTEGRATOR_WHFAST512 = 21,  // WHFast integrator, optimized for AVX512
        } integrator;
    enum {
        REB_BOUNDARY_NONE = 0,          // Do not check for anything (default)
        REB_BOUNDARY_OPEN = 1,          // Open boundary conditions. Removes particles if they leave the box 
        REB_BOUNDARY_PERIODIC = 2,      // Periodic boundary conditions
        REB_BOUNDARY_SHEAR = 3,         // Shear periodic boundary conditions, needs OMEGA variable
        } boundary;
    enum {
        REB_GRAVITY_NONE = 0,           // Do not calculate graviational forces
        REB_GRAVITY_BASIC = 1,          // Basic O(N^2) direct summation algorithm, choose this for shearing sheet and periodic boundary conditions
        REB_GRAVITY_COMPENSATED = 2,    // Direct summation algorithm O(N^2) but with compensated summation, slightly slower than BASIC but more accurate
        REB_GRAVITY_TREE = 3,           // Use the tree to calculate gravity, O(N log(N)), set opening_angle2 to adjust accuracy.
        REB_GRAVITY_MERCURIUS = 4,      // Special gravity routine only for MERCURIUS
        REB_GRAVITY_JACOBI = 5,         // Special gravity routine which includes the Jacobi terms for WH integrators 
        } gravity;

    // Datastructures for integrators
    struct reb_integrator_sei ri_sei;               // The SEI struct 
    struct reb_integrator_whfast ri_whfast;         // The WHFast struct 
    struct reb_integrator_whfast512 ri_whfast512;   // The WHFast512 struct 
    struct reb_integrator_saba ri_saba;             // The SABA struct 
    struct reb_integrator_ias15 ri_ias15;           // The IAS15 struct
    struct reb_integrator_mercurius ri_mercurius;   // The MERCURIUS struct
    struct reb_integrator_janus ri_janus;           // The JANUS struct 
    struct reb_integrator_eos ri_eos;               // The EOS struct 
    struct reb_integrator_bs ri_bs;                 // The BS struct

    // ODEs. Do not access these variables directly. Use functions provided instead.
    struct reb_ode** odes;      // all ode sets (includes nbody if BS is set as integrator)
    int N_odes;                 // number of ode sets
    int N_allocated_odes;
    int ode_warnings;

     // Callback functions
    void (*additional_forces) (struct reb_simulation* const r);             // Implement any additional (non-gravitational) forces here.
    void (*pre_timestep_modifications) (struct reb_simulation* const r);    // Executed just before eaach timestep. Used by REBOUNDx.
    void (*post_timestep_modifications) (struct reb_simulation* const r);   // Executed just after each timestep. Used by REBOUNDx.
    void (*heartbeat) (struct reb_simulation* r);                           // Executed at each timestep once. Use this to do extra output/work during a simulation.
    int (*key_callback) (struct reb_simulation* r, int key);                // Used when SERVER or OPENGL visualization is on. Gets called after completed timestep and if a key has been pressed. Return 1 if you want to skip default key commands.
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v);  // Allows for a velocity dependent coefficient of restitution (used for ring simulations)
    int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);        // Determines what happens when two particles collide.
    void (*free_particle_ap) (struct reb_particle* p);                      // Used by REBOUNDx.
    void (*extras_cleanup) (struct reb_simulation* r);                      // Used by REBOUNDx.
    void* extras;                                                           // Pointer to link to any additional (optional) libraries, e.g., REBOUNDx, ASSIST.
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// REBOUND API Functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Simulation life cycle

// Allocates memory for reb_simulation and initializes it.
DLLEXPORT struct reb_simulation* reb_simulation_create(void);
// Create a simulation object from a file. Set snapshot=-1 to load last snapshot.
DLLEXPORT struct reb_simulation* reb_simulation_create_from_file(char* filename, int64_t snapshot);
// Create a simulation object from a simulationarchive. Set snapshot=-1 to load last snapshot.
DLLEXPORT struct reb_simulation* reb_simulation_create_from_simulationarchive(struct reb_simulationarchive* sa, int64_t snapshot);
// Free simulation and all associated memory.
DLLEXPORT void reb_simulation_free(struct reb_simulation* const r);
// Only free memory in pointers of a simulation, but not the simulation itself.
DLLEXPORT void reb_simulation_free_pointers(struct reb_simulation* const r);
// Reset function pointers to default (NULL) values. Returns 1 if one ore more function pointers were not NULL before.
DLLEXPORT int reb_simulation_reset_function_pointers(struct reb_simulation* const r); 
// Reset all integrator variables.
DLLEXPORT void reb_simulation_reset_integrator(struct reb_simulation* r);
// Make a deep copy of simulation.
DLLEXPORT struct reb_simulation* reb_simulation_copy(struct reb_simulation* r);
// Compare r1 to r2. If exactly equal then 0 is returned, otherwise 1. If output_option=1, then difference is also printed on screen.
DLLEXPORT int reb_simulation_diff(struct reb_simulation* r1, struct reb_simulation* r2, int output_option);
// Setup simulation domain and root boxes. This needs to be called before particles are added if the tree code is used.
DLLEXPORT void reb_simulation_configure_box(struct reb_simulation* const r, const double boxsize, const int N_root_x, const int N_root_y, const int N_root_z); // Configure the boundary/root box

// Start webserver for visualization. Returns 0 on success.
DLLEXPORT int reb_simulation_start_server(struct reb_simulation* r, int port);
// Stop webserver.
DLLEXPORT void reb_simulation_stop_server(struct reb_simulation* r);

// Errors, warnings

// For fatal errors only. Print out an error message, then exit immediately and kill the process. Does no clean up memory.
DLLEXPORT void reb_exit(const char* const msg);
// Stop current integration in a nice way. Can be called from within heartbeat function.
DLLEXPORT void reb_simulation_stop(struct reb_simulation* const r);
// Print or store a warning message, then continue.
DLLEXPORT void reb_simulation_warning(struct reb_simulation* const r, const char* const msg);
// Print or store an error message, then continue.
DLLEXPORT void reb_simulation_error(struct reb_simulation* const r, const char* const msg);


// Output functions

// Write the simulation to file (simulationarchive format). Appends a snapshot if file exists.
DLLEXPORT void reb_simulation_save_to_file(struct reb_simulation* r, const char* filename);
// Schedule regular outputs to a file based on simulation time.
DLLEXPORT void reb_simulation_save_to_file_interval(struct reb_simulation* const r, const char* filename, double interval);
// Schedule regular outputs to a file based on wall time.
DLLEXPORT void reb_simulation_save_to_file_walltime(struct reb_simulation* const r, const char* filename, double walltime);
// Schedule regular outputs to a file based on number of steps taken.
DLLEXPORT void reb_simulation_save_to_file_step(struct reb_simulation* const r, const char* filename, uint64_t step);
// Write the simulation to a memory buffer (simulationarchive format).
DLLEXPORT void reb_simulation_save_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep);
// Output timing data to file. Appends file if it exists.
DLLEXPORT void reb_simulation_output_timing(struct reb_simulation* r, const double tmax);
// Output orbits to file. Appends file if it exists.
DLLEXPORT void reb_simulation_output_orbits(struct reb_simulation* r, char* filename);
// Output cartesian coordinates to file. Appends file if it exists.
DLLEXPORT void reb_simulation_output_ascii(struct reb_simulation* r, char* filename);
// Output velocity dispersion tensor to file. Used for ring simulations. Appends file if it exists.
DLLEXPORT void reb_simulation_output_velocity_dispersion(struct reb_simulation* r, char* filename);
// Function to allow for periodic outputs in heartbeat function. See examples on how to use it.
DLLEXPORT int reb_simulation_output_check(struct reb_simulation* r, double interval);


// Timestepping

// Advance simulation by 1 timestep.
DLLEXPORT void reb_simulation_step(struct reb_simulation* const r);
// Advance simulation by N_steps timesteps.
DLLEXPORT void reb_simulation_steps(struct reb_simulation* const r, unsigned int N_steps);
// Integrate simulation to at least time tmax (see exact_finish_time).
DLLEXPORT enum REB_STATUS reb_simulation_integrate(struct reb_simulation* const r, double tmax);
// Synchronize simulation if safe_mode is turned off by integrator to get physical coordinates.
DLLEXPORT void reb_simulation_synchronize(struct reb_simulation* r);


// Functions to operate on simulations

// Move the simulation to the heliocentric frame (particle with index 0 will be at the origin and at rest after calling this function).
DLLEXPORT void reb_simulation_move_to_hel(struct reb_simulation* const r);
// Move the simultion to the center of mass frame (the center of mass will be at the origin and at rest after calling this function).
DLLEXPORT void reb_simulation_move_to_com(struct reb_simulation* const r);
// Multiply x,y,z,vx,vy,vz of each particle in r with given scalars.
DLLEXPORT void reb_simulation_imul(struct reb_simulation* r, double scalar_pos, double scalar_vel);
// Add cartesian components of each particle of r2 to cartesian components in r. r2 and r must have same number of particles.
DLLEXPORT int reb_simulation_iadd(struct reb_simulation* r, struct reb_simulation* r2);
// Same as above but substract r2 from r component wise.
DLLEXPORT int reb_simulation_isub(struct reb_simulation* r, struct reb_simulation* r2);


// Diangnostic functions

// Return the sum of potential and kinetic energy
DLLEXPORT double reb_simulation_energy(const struct reb_simulation* const r);
// Returns the angular momentum.
DLLEXPORT struct reb_vec3d reb_simulation_angular_momentum(const struct reb_simulation* const r);
// Returns the center of mass of a simulation.
DLLEXPORT struct reb_particle reb_simulation_com(struct reb_simulation* r);
// Returns the center of mass of two particles.
DLLEXPORT struct reb_particle reb_particle_com_of_pair(struct reb_particle p1, struct reb_particle p2);
// Returns the center of mass of particles in the simulation within a given range.
DLLEXPORT struct reb_particle reb_simulation_com_range(struct reb_simulation* r, int first, int last);


// Functions to add and initialize particles

// Use this function to add particle pt to simulation r.
DLLEXPORT void reb_simulation_add(struct reb_simulation* const r, struct reb_particle pt);
// Use this function to add a particle to a simulation using orbital parameters or cartesian coordinates. See examples on usage.
DLLEXPORT void reb_simulation_add_fmt(struct reb_simulation* r, const char* fmt, ...);
// Same as reb_simulation_add_fmt() but returns the particle instead of adding it to the simulation. Still need simulation for G, center of mass, etc.
DLLEXPORT struct reb_particle reb_particle_from_fmt(struct reb_simulation* r, const char* fmt, ...);    
// This function sets up a Plummer sphere, N=number of particles, M=total mass, R=characteristic radius. Particles get added to simulation.
DLLEXPORT void reb_simulation_add_plummer(struct reb_simulation* r, int _N, double M, double R); 
// Returns a particle given a set of orbital parameters. Also sets err to error code if initialization failed.
DLLEXPORT struct reb_particle reb_particle_from_orbit_err(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f, int* err);
// Same as above but without error code.
DLLEXPORT struct reb_particle reb_particle_from_orbit(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f);
// Returns a particle given a set of Pal orbital parameters.
DLLEXPORT struct reb_particle reb_particle_from_pal(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy);
// Returns a reb_particle structure with fields/hash/ptrs initialized to nan/0/NULL. 
DLLEXPORT struct reb_particle reb_particle_nan(void);


// Functions to access, remove, and operate on particles

// Remove all particles
DLLEXPORT void reb_simulation_remove_all_particles(struct reb_simulation* const r);
// Remove one particle. keep_sorted flag can be set to 1 to maintain order of remaining particles.
DLLEXPORT int reb_simulation_remove_particle(struct reb_simulation* const r, int index, int keep_sorted);
// Remove one particle. Use hash to find particle.
DLLEXPORT int reb_simulation_remove_particle_by_hash(struct reb_simulation* const r, uint32_t hash, int keep_sorted);
// Returns pointer to particle with given hash.
DLLEXPORT struct reb_particle* reb_simulation_particle_by_hash(struct reb_simulation* const r, uint32_t hash);
// Same as above but searches on all nodes if MPI is enabled. Returns copy instead of a pointer because particle might be on a different node.
DLLEXPORT struct reb_particle reb_simulation_particle_by_hash_mpi(struct reb_simulation* const r, uint32_t hash);
// Returns a particle's index in the simulation given a pointer to the particle. Returns -1 if not found. 
DLLEXPORT int reb_simulation_particle_index(struct reb_particle* p); 
// Subtract particle p2 from p1
DLLEXPORT void reb_particle_isub(struct reb_particle* p1, struct reb_particle* p2);
// Add particle p2 to p1
DLLEXPORT void reb_particle_iadd(struct reb_particle* p1, struct reb_particle* p2);
// Multiply x,y,z,vx,vy,vz,m of p1 with given value.
DLLEXPORT void reb_particle_imul(struct reb_particle* p1, double value);
// Return the distance between particle p1 and p2
DLLEXPORT double reb_particle_distance(struct reb_particle* p1, struct reb_particle* p2);
// Compares two particles, ignoring pointers. Returns 1 if particles differ, 0 if they are exactly equal.
DLLEXPORT int reb_particle_diff(struct reb_particle p1, struct reb_particle p2); 


// Chaos indicators

// Turn on MEGNO/Lyapunov calculation. Uses random seen in simulation.
DLLEXPORT void reb_simulation_init_megno(struct reb_simulation* const r);
// Same as above but used given random seend. Useful to reproduce same results every time.
DLLEXPORT void reb_simulation_init_megno_seed(struct reb_simulation* const r, unsigned int seed);
// Returns the current MEGNO value,
DLLEXPORT double reb_simulation_megno(struct reb_simulation* r);
// Returns the largest Lyapunov characteristic number (LCN).
DLLEXPORT double reb_simulation_lyapunov(struct reb_simulation* r);


// Built in mercurius switching functions

DLLEXPORT double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit);  // default
DLLEXPORT double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit);
DLLEXPORT double reb_integrator_mercurius_L_C4(const struct reb_simulation* const r, double d, double dcrit);
DLLEXPORT double reb_integrator_mercurius_L_C5(const struct reb_simulation* const r, double d, double dcrit);


// Built in collision resolve functions

DLLEXPORT int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c); // halts a simulation when a collision occurs
DLLEXPORT int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);
DLLEXPORT int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);


// Random sampling - These functions only use the simulation object for a seed. If r=NULL time and PID are used as a seed.

DLLEXPORT double reb_random_uniform(struct reb_simulation* r, double min, double max);
DLLEXPORT double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
DLLEXPORT double reb_random_normal(struct reb_simulation* r, double variance);
DLLEXPORT double reb_random_rayleigh(struct reb_simulation* r, double sigma);


// Miscellaneous functions

// Calculate a hash value for a string.
DLLEXPORT uint32_t reb_hash(const char* str);   
// Returns the angle f wrapped in the interval from 0 to 2*pi
DLLEXPORT double reb_mod2pi(double f);
// True anomaly for a given eccentricity and mean anomaly
DLLEXPORT double reb_M_to_f(double e, double M);
// True anomaly for a given eccentricity and eccentric anomaly
DLLEXPORT double reb_E_to_f(double e, double M);
// Eccentric anomaly for a given eccentricity and mean anomaly
DLLEXPORT double reb_M_to_E(double e, double M);


// Simulationarchive

// Simulationarchive structure
struct reb_simulationarchive{
    FILE* inf;                      // File pointer (will be kept open)
    char* filename;                 // Filename of open file. This is NULL if this is a memory-mapped file (using fmemopen)
    int version;                    // Simulationarchive version
    double auto_interval;           // Interval setting used to create SA (if used)
    double auto_walltime;           // Walltime setting used to create SA (if used)
    uint64_t auto_step;   // Steps in-between SA snapshots (if used)
    int64_t nblobs;                    // Total number of snapshots (including initial binary)
    uint64_t* offset;               // Index of offsets in file (length nblobs)
    double* t;                      // Index of simulation times in file (length nblobs)
};
// Allocate memory for a simulationarchive and initialize it with a file.
DLLEXPORT struct reb_simulationarchive* reb_simulationarchive_create_from_file(const char* filename);
// Free memory allocated by simulationarchive
DLLEXPORT void reb_simulationarchive_free(struct reb_simulationarchive* sa);


// Orbit calculation

// Structure containing orbital elements for Keplerian orbits.
struct reb_orbit {
    double d;               // Radial distance from central object
    double v;               // velocity relative to central object's velocity
    double h;               // Specific angular momentum
    double P;               // Orbital period
    double n;               // Mean motion
    double a;               // Semi-major axis
    double e;               // Eccentricity
    double inc;             // Inclination
    double Omega;           // Longitude of ascending node
    double omega;           // Argument of pericenter
    double pomega;          // Longitude of pericenter
    double f;               // True anomaly
    double M;               // Mean anomaly
    double l;               // Mean Longitude
    double theta;           // True Longitude
    double T;               // Time of pericenter passage
    double rhill;           // Circular Hill radius 
    double pal_h;           // Cartesian component of the eccentricity, h = e*sin(pomega) 
    double pal_k;           // Cartesian component of the eccentricity, k = e*cos(pomega) 
    double pal_ix;          // Cartesian component of the inclination, ix = 2*sin(i/2)*cos(Omega)
    double pal_iy;          // Cartesian component of the inclination, ix = 2*sin(i/2)*sin(Omega)
    struct reb_vec3d hvec;  // specific angular momentum vector
    struct reb_vec3d evec;  // eccentricity vector (mag=ecc, points toward peri)
};
// Calculates all orbital elements of the particle p, assuming gravitational constant G and the given primary.
DLLEXPORT struct reb_orbit reb_orbit_from_particle(double G, struct reb_particle p, struct reb_particle primary);


// ODE functions

// Defines one Ordinary Differential Equation (ODE) so that it can be integrated with REBOUND
struct reb_ode{
    unsigned int length;        // number of components / dimenion
    double* y;                  // Pointer to current state 
    unsigned int needs_nbody;   // 1: ODE needs N-body particles to calculate RHS
    void* ref;                  // Optional pointer to any additional data needed for derivative calculation
    void (*derivatives)(struct reb_ode* const ode, double* const yDot, const double* const y, const double t);  // Function pointer to right hand side of ODE
    void (*getscale)(struct reb_ode* const ode, const double* const y0, const double* const y1);                // Function pointer, sets scales for components (optional) 
    void (*pre_timestep)(struct reb_ode* const ode, const double* const y0);                                    // Function pointer, gets called just before the ODE integration (optional)
    void (*post_timestep)(struct reb_ode* const ode, const double* const y0);                                   // Function pointer, gets called just after the ODE integration (optional)

    // Internal use
    unsigned int N_allocated;   
    double* scale;
    double* C;                  // Temporary internal array (extrapolation) 
    double** D;                 // Temporary internal array (extrapolation) 
    double* y1;                 // Temporary internal array (state during the step) 
    double* y0Dot;              // Temporary internal array (derivatives at beginning of step)
    double* yDot;               // Temporary internal array (derivatives)
    double* yTmp;               // Temporary internal array (midpoint method)
    struct reb_simulation* r;   // weak reference to main simulation 
};
// Allocate memory for an ODE struct, initialize it, and attach it to the simulation. See examples for details on usage.
DLLEXPORT struct reb_ode* reb_ode_create(struct reb_simulation* r, unsigned int length);
// Free an ODE struct.
DLLEXPORT void reb_ode_free(struct reb_ode* ode);


// Variational equations

// Struct describing the properties of a set of variational equations.
// If testparticle is set to -1, then it is assumed that all particles are massive
// and all particles influence all other particles. If testparticle is >=0 then 
// the particle with that index is assumed to be a testparticle, i.e. it does not 
// influence other particles. For second order variational equation, index_1st_order_a/b 
// is the index in the particle array that corresponds to the 1st order variational 
// equations.
struct reb_variational_configuration{
    struct reb_simulation* sim; // Reference to the simulation.
    int order;                  // Order of the variational equation. 1 or 2. 
    int index;                  // Index of the first variational particle in the particles array.
    int testparticle;           // Is this variational configuration describe a test particle? -1 if not.
    int index_1st_order_a;      // Used for 2nd order variational particles only: Index of the first order variational particle in the particles array.
    int index_1st_order_b;      // Used for 2nd order variational particles only: Index of the first order variational particle in the particles array.
    double lrescale;            // Accumulates the logarithm of rescalings
};

// Add and initialize a set of first order variational particles
// If testparticle is >= 0, then only one variational particle (the test particle) will be added.
// If testparticle is -1, one variational particle for each real particle will be added.
// Returns the index of the first variational particle added
DLLEXPORT int reb_simulation_add_variation_1st_order(struct reb_simulation* const r, int testparticle);

// Add and initialize a set of second order variational particles
// Note: a set of second order variational particles requires two sets of first order variational equations.
// If testparticle is >= 0, then only one variational particle (the test particle) will be added.
// If testparticle is -1, one variational particle for each real particle will be added.
// index_1st_order_a is the index of the corresponding first variational particles.
// index_1st_order_b is the index of the corresponding first variational particles.
// Returns the index of the first variational particle added
DLLEXPORT int reb_simulation_add_variation_2nd_order(struct reb_simulation* const r, int testparticle, int index_1st_order_a, int index_1st_order_b);

// Rescale all sets of variational particles if their size gets too large (>1e100).
// This can prevent an overflow in floating point numbers. The logarithm of the rescaling
// factor is stored in the reb_variational_configuration's lrescale variable. 
// This function is called automatically every timestep. To avoid automatic rescaling,
// set the reb_variational_configuration's lrescale variable to -1.
// For this function to work, the positions and velocities needs to be synchronized. 
// A warning is presented if the integrator is not synchronized. 
DLLEXPORT void reb_simulation_rescale_var(struct reb_simulation* const r);

// These functions calculates the first/second derivative of a Keplerian orbit. 
//   Derivatives of Keplerian orbits are required for variational equations, in particular
//   for optimization problems. 
//   The derivative is calculated with respect to the variables that appear in the function name.
//   One variable implies that a first derivative is returned, two variables implies that a second
//   derivate is returned. Classical orbital parameters and those introduced by Pal (2009) are 
//   supported. Pal coordinates have the advantage of being analytical (i.e. infinite differentiable).
//   Classical orbital parameters may have singularities, for example when e is close to 0.
//   Note that derivatives with respect to Cartesian coordinates are trivial and therefore not
//   implemented as seperate functions. 
//   The following variables are supported: a, e, inc, f, omega, Omega, h, k, ix, iy and m (mass). 
// The functions return the derivative as a particle structre. Each structure element is a derivative.
// The paramter po is the original particle for which the derivative is to be calculated.
DLLEXPORT struct reb_particle reb_particle_derivative_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_h(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k_k(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_h_h(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_lambda_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_h_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k_h(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_a(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_ix_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_iy_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_h_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_lambda_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_lambda_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_h_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_k_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_ix_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_h(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_k(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_a(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_lambda(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_h(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_k(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_ix(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_iy(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_m(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e_e(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_inc(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_inc_inc(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_Omega_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_omega_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_f_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_e(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_inc(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_a_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e_inc(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_e_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_e(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_inc_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_inc_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_inc_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_inc(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_omega_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_Omega_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_Omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_omega_f(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_omega(double G, struct reb_particle primary, struct reb_particle po);
DLLEXPORT struct reb_particle reb_particle_derivative_m_f(double G, struct reb_particle primary, struct reb_particle po);


// Functions to convert between coordinate systems

// Rotations
DLLEXPORT struct reb_rotation reb_rotation_inverse(const struct reb_rotation q);
DLLEXPORT struct reb_rotation reb_rotation_mul(const struct reb_rotation p, const struct reb_rotation q);

DLLEXPORT struct reb_rotation reb_rotation_identity();
DLLEXPORT struct reb_rotation reb_rotation_normalize(const struct reb_rotation q);
DLLEXPORT struct reb_rotation reb_rotation_conjugate(const struct reb_rotation q);
DLLEXPORT struct reb_rotation reb_rotation_init_angle_axis(const double angle, struct reb_vec3d axis);
DLLEXPORT struct reb_rotation reb_rotation_init_from_to(struct reb_vec3d from, struct reb_vec3d to);
DLLEXPORT struct reb_rotation reb_rotation_init_orbit(const double Omega, const double inc, const double omega);
DLLEXPORT struct reb_rotation reb_rotation_init_to_new_axes(struct reb_vec3d newz, struct reb_vec3d newx);
DLLEXPORT struct reb_rotation reb_rotation_slerp(struct reb_rotation q1, struct reb_rotation q2, double t);

// transformations to/from vec3d
DLLEXPORT struct reb_vec3d reb_tools_spherical_to_xyz(const double mag, const double theta, const double phi);
DLLEXPORT void reb_tools_xyz_to_spherical(struct reb_vec3d const xyz, double* mag, double* theta, double* phi);

DLLEXPORT struct reb_vec3d reb_vec3d_mul(const struct reb_vec3d v, const double s);
DLLEXPORT struct reb_vec3d reb_vec3d_add(const struct reb_vec3d v, const struct reb_vec3d w);
DLLEXPORT double reb_vec3d_length_squared(const struct reb_vec3d v);
DLLEXPORT double reb_vec3d_dot(const struct reb_vec3d a, const struct reb_vec3d b);
DLLEXPORT struct reb_vec3d reb_vec3d_cross(const struct reb_vec3d a, const struct reb_vec3d b);
DLLEXPORT struct reb_vec3d reb_vec3d_normalize(const struct reb_vec3d v);
DLLEXPORT struct reb_vec3d reb_vec3d_rotate(struct reb_vec3d v, const struct reb_rotation q);
DLLEXPORT void reb_vec3d_irotate(struct reb_vec3d *v, const struct reb_rotation q);
DLLEXPORT void reb_particle_irotate(struct reb_particle* p, const struct reb_rotation q);
DLLEXPORT void reb_simulation_irotate(struct reb_simulation* const sim, const struct reb_rotation q);

DLLEXPORT void reb_rotation_to_orbital(struct reb_rotation q, double* Omega, double* inc, double* omega);

#ifdef MPI
void reb_mpi_init(struct reb_simulation* const r);
void reb_mpi_finalize(struct reb_simulation* const r);
#endif // MPI

#ifdef OPENMP
// Wrapper method to set number of OpenMP threads from python.
DLLEXPORT void reb_omp_set_num_threads(int num_threads);
#endif // OPENMP

// The following stuctures are related to OpenGL/WebGL visualization. Nothing to be changed by the user.

struct reb_orbit_opengl {
    float x,y,z;
    float a, e, f;
    float omega, Omega, inc;
};

struct reb_vec3df {
    float x,y,z;
};

struct reb_vec4df {
    float x,y,z,r;
};

struct reb_server_data {
    struct reb_simulation* r;
    struct reb_simulation* r_copy;
    int port;
    int need_copy;
    int ready;
#ifdef SERVER
#ifdef _WIN32
    SOCKET socket;
    HANDLE mutex;          // Mutex to allow for copying
#else // _WIN32
    int socket;
    pthread_mutex_t mutex;          // Mutex to allow for copying
    pthread_t server_thread;
#endif // _WIN32
#endif // SERVER
};

struct reb_display_settings {
    struct reb_mat4df view;
    int spheres;                    // Switches between point sprite and real spheres.
    int pause;                      // Pauses visualization, but keep simulation running
    int wire;                       // Shows/hides orbit wires.
    int past;                       // Shows/hides past particle positions.
    unsigned int past_N;            // Number of past particle positions.
    int onscreentext;               // Shows/hides onscreen text.
    int onscreenhelp;               // Shows/hides onscreen help.
    int multisample;                // Turn off/on multisampling.
    int ghostboxes;                 // Shows/hides ghost boxes.
    int reference;                  // reb_particle used as a reference for centering.
};

struct reb_display_data {
    struct reb_display_settings s;
    struct reb_simulation* r;
    struct reb_simulation* r_copy;
    struct reb_vec4df* particle_data;
    struct reb_orbit_opengl* orbit_data;
    uint64_t N_allocated;
    double mouse_x;
    double mouse_y;
    double retina;
#ifndef _WIN32
    int need_copy;
    pthread_mutex_t mutex;          // Mutex to allow for copying
    pthread_t compute_thread;
#endif // _WIN32
#ifdef __EMSCRIPTEN__
    int connection_status;
#endif
    uint64_t past_last_steps_done;
    unsigned int past_N_allocated;
    unsigned int past_current_index;
    unsigned int mouse_action;      
    unsigned int key_mods;      
    unsigned int particle_buffer;
    unsigned int particle_buffer_current;
    unsigned int orbit_buffer;
    unsigned int orbit_buffer_current;
    void* window;
    struct {
        unsigned int texture;
        unsigned int program;
        unsigned int vao;
        unsigned int pos_location;
        unsigned int ypos_location;
        unsigned int scale_location;
        unsigned int aspect_location;
        unsigned int screen_aspect_location;
        unsigned int charval_buffer;
    } shader_simplefont;
    struct {
        unsigned int program;
        unsigned int box_vao;
        unsigned int cross_vao;
        unsigned int mvp_location;
        unsigned int color_location;
    } shader_box;
    struct {
        unsigned int mvp_location;
        unsigned int color_location;
        unsigned int current_vertex_location;
        unsigned int past_N_location;
        unsigned int N_real_location;
        unsigned int program;
        unsigned int particle_vao;
    } shader_point;
    struct {
        unsigned int mvp_location;
        unsigned int program;
        unsigned int particle_vao_current;
        unsigned int particle_vao;
    } shader_sphere;
    struct {
        unsigned int mvp_location;
        unsigned int current_vertex_location;
        unsigned int past_N_location;
        unsigned int N_real_location;
        unsigned int program;
        unsigned int particle_vao_current;
        unsigned int particle_vao;
        unsigned int vertex_count;
    } shader_orbit;

};

// Display settings initialization
DLLEXPORT void reb_simulation_add_display_settings(struct reb_simulation* r);

// Matrix methods
DLLEXPORT struct reb_mat4df reb_mat4df_identity();
DLLEXPORT struct reb_mat4df reb_mat4df_scale(struct reb_mat4df m, float x, float y, float z);
DLLEXPORT void reb_mat4df_print(struct reb_mat4df m);
DLLEXPORT int reb_mat4df_eq(struct reb_mat4df A, struct reb_mat4df B);
DLLEXPORT struct reb_vec3df reb_mat4df_get_scale(struct reb_mat4df m);
DLLEXPORT struct reb_mat4df reb_mat4df_translate(struct reb_mat4df m, float x, float y, float z);
DLLEXPORT struct reb_mat4df reb_mat4df_multiply(struct reb_mat4df A, struct reb_mat4df B);
DLLEXPORT struct reb_mat4df reb_rotation_to_mat4df(struct reb_rotation A);
DLLEXPORT struct reb_mat4df reb_mat4df_ortho(float l, float r, float b, float t, float n, float f);


// Declarations and functions needed internally or by python interface only.
void reb_sigint_handler(int signum);

// Used in the binary file to identify data blobs
struct reb_simulationarchive_blob {  
    int32_t index;                   // Index of previous blob (binary file is 0, first blob is 1)
    int32_t offset_prev;             // Offset to beginning of previous blob (size of previous blob).
    int32_t offset_next;             // Offset to end of following blob (size of following blob).
};
// Binary field descriptors are used to identify data blobs in simulationarchives.
struct reb_binary_field_descriptor {
    uint32_t type;          // Unique id for each field. Should not change between versions. Ids should not be reused.
    enum {
        REB_DOUBLE = 0,
        REB_INT = 1,
        REB_UINT = 2,                // Same as UINT32
        REB_UINT32 = 3,
        REB_INT64 = 4,
        REB_UINT64 = 5,
        // REB_ULONGLONG = 6,        // No longer used. Using explicit lengths instead.
        REB_VEC3D = 7,
        REB_PARTICLE = 8,
        REB_POINTER = 9,
        REB_POINTER_ALIGNED = 10,    // memory aligned to 64 bit boundary for AVX512
        REB_DP7 = 11,                // Special datatype for IAS15
        REB_OTHER = 12,              // Fields that need special treatment during input and/or output
        REB_FIELD_END = 13,          // Special type to indicate end of blob
        REB_FIELD_NOT_FOUND = 14,    // Special type used to throw error messages
        REB_PARTICLE4 = 15,          // Used for WHFast512
        REB_POINTER_FIXED_SIZE = 16, // A pointer with a fixed size.
    } dtype;
    char name[1024];
    size_t offset;              // Offset of the storage location relative to the beginning of reb_simulation
    size_t offset_N;            // Offset of the storage location for the size relative to the beginning of reb_simulation
    size_t element_size;        // Size in bytes of each element (only used for pointers, dp7, etc)
};
DLLEXPORT extern const struct reb_binary_field_descriptor reb_binary_field_descriptor_list[]; // List of blobs. Implemented in output.c
DLLEXPORT struct reb_binary_field_descriptor reb_binary_field_descriptor_for_type(int type);
DLLEXPORT struct reb_binary_field_descriptor reb_binary_field_descriptor_for_name(const char* name);

// Possible errors that might occur during binary file reading.
enum reb_simulation_binary_error_codes {
    REB_SIMULATION_BINARY_WARNING_NONE = 0,
    REB_SIMULATION_BINARY_ERROR_NOFILE = 1,
    REB_SIMULATION_BINARY_WARNING_VERSION = 2,
    REB_SIMULATION_BINARY_WARNING_POINTERS = 4,
    REB_SIMULATION_BINARY_WARNING_PARTICLES = 8,
    REB_SIMULATION_BINARY_ERROR_FILENOTOPEN = 16,
    REB_SIMULATION_BINARY_ERROR_OUTOFRANGE = 32,
    REB_SIMULATION_BINARY_ERROR_SEEK = 64,
    REB_SIMULATION_BINARY_WARNING_FIELD_UNKOWN = 128,
    REB_SIMULATION_BINARY_ERROR_INTEGRATOR = 256,
    REB_SIMULATION_BINARY_WARNING_CORRUPTFILE = 512,
    REB_SIMULATION_BINARY_ERROR_OLD = 1024,
};


struct reb_binary_field { // This structure is used to save and load binary files.
    uint32_t type;  // type as given by reb_binary_field_descriptor
    uint64_t size;  // Size in bytes of field (only counting what follows, not the binary field, itself).
};

DLLEXPORT void reb_simulation_init(struct reb_simulation* r); // Used internally and by python. Should not be called by the user.
DLLEXPORT void reb_simulation_update_acceleration(struct reb_simulation* r); // Used by REBOUNDx
DLLEXPORT void reb_simulation_update_tree(struct reb_simulation* const r);
DLLEXPORT int reb_simulation_get_next_message(struct reb_simulation* const r, char* const buf); // Get the next stored warning message. Used only if save_messages==1. Return value is 0 if no messages are present, 1 otherwise.
DLLEXPORT int reb_check_fp_contract(); // Returns 1 if floating point contraction are enabled. 0 otherwise.
DLLEXPORT size_t reb_simulation_struct_size();
DLLEXPORT char* reb_simulation_diff_char(struct reb_simulation* r1, struct reb_simulation* r2); // Return the difference between two simulations as a human readable difference. Returned pointer needs to be freed.
DLLEXPORT void reb_simulation_set_collision_resolve(struct reb_simulation* r, int (*resolve) (struct reb_simulation* const r, struct reb_collision c)); // Used from python 
DLLEXPORT void reb_simulation_get_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // NULL pointers will not be set.
DLLEXPORT void reb_simulation_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // Null pointers will be ignored.
DLLEXPORT void reb_simulation_output_free_stream(char* buf);
DLLEXPORT struct reb_particle reb_simulation_jacobi_com(struct reb_particle* p); // Returns the Jacobi center of mass for a given particle. Used by python. Particle needs to be in a simulation.
DLLEXPORT struct reb_orbit reb_orbit_from_particle_err(double G, struct reb_particle p, struct reb_particle primary, int* err);
DLLEXPORT void reb_simulation_create_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, int64_t snapshot, enum reb_simulation_binary_error_codes* warnings);
DLLEXPORT void reb_simulation_copy_with_messages(struct reb_simulation* r_copy,  struct reb_simulation* r, enum reb_simulation_binary_error_codes* warnings); // used from python
DLLEXPORT void reb_simulationarchive_init_from_buffer_with_messages(struct reb_simulationarchive* sa, char* buf, size_t size, struct reb_simulationarchive* sa_index, enum reb_simulation_binary_error_codes* warnings);
DLLEXPORT void reb_simulationarchive_create_from_file_with_messages(struct reb_simulationarchive* sa, const char* filename,  struct reb_simulationarchive* sa_index, enum reb_simulation_binary_error_codes* warnings);
DLLEXPORT void reb_simulationarchive_free_pointers(struct reb_simulationarchive* sa);

DLLEXPORT void reb_particles_transform_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active); // p_mass: Should be the same particles array as ps for real particles. If passing variational particles in ps, p_mass should be the corresponding array of real particles.
DLLEXPORT void reb_particles_transform_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j,const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const unsigned int N_active);

// Democratic heliocentric coordinates
DLLEXPORT void reb_particles_transform_inertial_to_democraticheliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_democraticheliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_democraticheliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);

// WHDS
DLLEXPORT void reb_particles_transform_inertial_to_whds_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_whds_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);
DLLEXPORT void reb_particles_transform_whds_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const unsigned int N_active);


// Temporary. Function declarations needed by REBOUNDx 
DLLEXPORT void reb_integrator_ias15_reset(struct reb_simulation* r);         ///< Internal function used to call a specific integrator
DLLEXPORT void reb_integrator_ias15_part2(struct reb_simulation* r);         ///< Internal function used to call a specific integrator
DLLEXPORT void reb_integrator_whfast_from_inertial(struct reb_simulation* const r);   ///< Internal function to the appropriate WHFast coordinates from inertial
DLLEXPORT void reb_integrator_whfast_to_inertial(struct reb_simulation* const r); ///< Internal function to move back from particular WHFast coordinates to inertial
DLLEXPORT void reb_integrator_whfast_reset(struct reb_simulation* r);		///< Internal function used to call a specific integrator
DLLEXPORT int reb_integrator_whfast_init(struct reb_simulation* const r);    ///< Internal function to check errors and allocate memory if needed
DLLEXPORT void reb_whfast_interaction_step(struct reb_simulation* const r, const double _dt);///< Internal function
DLLEXPORT void reb_whfast_jump_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
DLLEXPORT void reb_whfast_kepler_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
DLLEXPORT void reb_whfast_com_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
#endif // _MAIN_H
