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

// Operating system specific options.
// Windows requires special treatment.
#ifdef _WIN64
#define _LP64   // automatically defined on 64bit Linux and MacOS
#endif // _WIN64
#ifdef _WIN32
#define _USE_MATH_DEFINES // Windows (MVSC) does not include math constants by default.
#define _NO_CRT_STDIO_INLINE // WIN32 to use _vsprintf_s
#include <WinSock2.h>
#define _WINSOCKAPI_ //stops windows.h including winsock.h
#include <windows.h>
#define REB_RESTRICT
#define DLLEXPORT __declspec(dllexport)
#define __restrict__
#ifdef _MSC_VER
#pragma comment(lib, "legacy_stdio_definitions.lib") // for printf, etc
#endif
#elif __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#define REB_RESTRICT
#define DLLEXPORT
#else // Linux and MacOS
#define REB_RESTRICT restrict
#define DLLEXPORT
#endif // _WIN32

#include <stdlib.h> // for size_t
#include <stdint.h> // for integer types


// Forward declarations
struct reb_simulationarchive;   // Opaque pointer. Implemented in simulationarchive.h
struct reb_server_data;         // Opaque pointer. Implemented in server.h
struct reb_treecell;            // Opaque pointer. Implemented in tree.h
struct reb_display_data;        // Opaque pointer. Implemented in display.h
struct reb_particle_int;        // Opaque pointer. Implemented in integrator_janus.h
struct reb_name_hash_item;      // Opaque pointer. Implemented in particle.h
struct reb_simulation;          // Implemented below.
struct reb_display_settings;    // Implemented below.
struct reb_variational_configuration;  // Implemented below.

// Particle structure
// Note: Size is 128 bytes, which corresponds to two L1 cache lines (one on Apple silicon).
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
    double m;                   // Mass in code units
    double r;                   // Physical radius in code units
    const char* name;           // Pointer to a NULL terminated string with the particle's name.
#if !defined(_LP64)
    char pad2[4];               // Padding. ap is not padded to 8 bytes
#endif
    void* ap;                   // This pointer allows REBOUNDx to add additional properties to the particle.
#if !defined(_LP64)
    char pad3[4];               // Padding. ap is short by 4 bytes
#endif
    struct reb_simulation* sim; // Pointer to the parent simulation.
#if !defined(_LP64)
    char pad4[4];               // Padding. sim is short by 4 bytes
#endif
};

// Generic 3d vector
struct reb_vec3d {
    double x;
    double y;
    double z;
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

// Structure representing one particle-particle collision
struct reb_collision{
    size_t p1;              // Index of first particle involved in collision
    size_t p2;              // Index of second particle
    struct reb_vec6d gb;    // Offset due to boundary conditions
    size_t ri;              // Root cell index (MPI only)
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
    double epsilon;                 // Precision control parameter
    double min_dt;                  // Minimal timestep
    enum {
        REB_IAS15_INDIVIDUAL = 0,   // fractional error is calculated separately for each particle
        REB_IAS15_GLOBAL = 1,       // fractional error is calculated globally (was default until 01/2024)
        REB_IAS15_PRS23 = 2,        // Pham, Rein & Spiegel (2023) timestep criterion (default since 01/2024)
        REB_IAS15_AARSETH85 = 3,    // Aarseth (1985) timestep criterion
    } adaptive_mode;
    uint64_t iterations_max_exceeded; // Counter how many times the iteration did not converge. 
    size_t N_allocated;          
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
    enum {
        REB_MERCURIUS_MODE_WH = 0,
        REB_MERCURIUS_MODE_ENCOUNTER = 1,
    } mode;
    size_t encounter_N;             // Number of particles currently having an encounter
    size_t encounter_N_active;      // Number of active particles currently having an encounter
    unsigned int tponly_encounter;  // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles
    size_t N_allocated;
    size_t N_allocated_additional_forces;
    size_t N_allocated_dcrit;       // Current size of dcrit arrays
    double* dcrit;                  // Precalculated switching radii for particles
    struct reb_particle* REB_RESTRICT particles_backup; //  contains coordinates before Kepler step for encounter prediction
    struct reb_particle* REB_RESTRICT particles_backup_additional_forces; // contains coordinates before Kepler step for encounter prediction
    size_t* encounter_map;          // Map to represent which particles are integrated with ias15
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

// Leapfrog Integrator (TU splitting)
struct reb_integrator_leapfrog {
    unsigned int order;
};

// TRACE (Lu Hernandez & Rein 2024)
struct reb_integrator_trace {
    int (*S) (struct reb_simulation* const r, const size_t i, const size_t j);
    int (*S_peri) (struct reb_simulation* const r, const size_t j);

    enum {
        REB_TRACE_PERI_PARTIAL_BS = 0,
        REB_TRACE_PERI_FULL_BS = 1,
        REB_TRACE_PERI_FULL_IAS15 = 2,
    } peri_mode;

    double r_crit_hill;
    double peri_crit_eta;

    // Internal use
    enum {
        REB_TRACE_MODE_INTERACTION = 0, // Interaction step
        REB_TRACE_MODE_KEPLER = 1,      // Kepler step
        REB_TRACE_MODE_NONE = 2,        // In-between steps, to avoid calculate_accelerations
        REB_TRACE_MODE_FULL = 3,        // Doing everything in one step
    } mode;
    size_t encounter_N;                 // Number of particles currently having an encounter
    size_t encounter_N_active;          // Number of active particles currently having an encounter

    size_t N_allocated;
    size_t N_allocated_additional_forces;
    unsigned int tponly_encounter;      // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles

    struct reb_particle* REB_RESTRICT particles_backup; //  Contains coordinates before the entire step
    struct reb_particle* REB_RESTRICT particles_backup_kepler; //  Contains coordinates before kepler step
    struct reb_particle* REB_RESTRICT particles_backup_additional_forces; // For additional forces

    size_t* encounter_map;              // Map to represent which particles are integrated with BS
    size_t* encounter_map_backup;       // Contains encounter map from after pre-ts check. Used to retain memory of CEs flagged at this step.
    struct reb_vec3d com_pos;           // Used to keep track of the centre of mass during the timestep
    struct reb_vec3d com_vel;

    int* current_Ks; // Tracking K_ij for the entire timestep
    unsigned int current_C; // Tracking C for the entire timestep
    unsigned int force_accept; // Force accept for irreversible steps: collisions and adding particles
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
        REB_WHFAST_COORDINATES_BARYCENTRIC = 3,                 // Barycentric coordinates
    } coordinates;                                              // Coordinate system used in Hamiltonian splitting
    unsigned int recalculate_coordinates_this_timestep;         // 1: recalculate coordinates from inertial coordinates
    unsigned int safe_mode;                                     // 0: Drift Kick Drift scheme (default), 1: combine first and last sub-step.
    unsigned int keep_unsynchronized;                           // 1: continue from unsynchronized state after synchronization 

    // Internal use
    size_t N_allocated;
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    size_t N_allocated_var;
    struct reb_particle* REB_RESTRICT p_jh_var; // Jacobi coordinates for variational equations
    size_t N_allocated_temp;
    struct reb_particle* REB_RESTRICT p_temp;   // Used for lazy implementer's kernel 
    unsigned int is_synchronized;
    unsigned int timestep_warning;
    unsigned int recalculate_coordinates_but_not_synchronized_warning;
};

// Special particle struct for WHFast512
struct reb_particle_avx512; // Opaque pointer. Implemented in integrator_whfast.h

// WHFast512 Integrator (Javaheri & Rein 2023)
struct reb_integrator_whfast512 {
    unsigned int gr_potential;          // 1: Turn on GR potential of central object, 0 (default): no GR potential
    unsigned int N_systems;             // Number of systems to be integrator in parallel: 1 (default, up to 8 planets), 2 (up to 4 planets each), 4 (2 planets each)
    unsigned int keep_unsynchronized;   // 1: continue from unsynchronized state after synchronization 

    // Internal use
    unsigned int is_synchronized;
    size_t N_allocated;
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

// Available return values for collision resolve functions
enum REB_COLLISION_RESOLVE_OUTCOME {
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_NONE = 0,
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P1 = 1,
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_P2 = 2,
    REB_COLLISION_RESOLVE_OUTCOME_REMOVE_BOTH = 3,
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

// Janus integrator (Rein & Tamayo 2018)
struct reb_integrator_janus {
    double scale_pos;       // Scale of position grid. Default 1e-16
    double scale_vel;       // Scale of velocity grid. Default 1e-16
    unsigned int order;     // Order: 2 (default), 4, 6, 8, 10 
    unsigned int recalculate_integer_coordinates_this_timestep;  // Set to 1 if particles have been modified

    // Internal use
    struct reb_particle_int* REB_RESTRICT p_int;
    size_t N_allocated;
};

// Generic custom integrator
struct reb_integrator {
    void (*step)(struct reb_simulation* r);         // Performs one timestep. Timestep should be r->dt for a non-adaptive integrator. Need to update r->t in this routine.
    void (*synchronize)(struct reb_simulation* r);  // Synchronizes particle state. Optional. Set to NULL if not used.
    void (*reset)(struct reb_simulation* r);        // Reset intergrator state to default and free all memory. Optional. Set to NULL if not used.
    void* data;                                     // Pointer to any internal data/memory required by the integrator. Optional.
    size_t data_size;                               // Size of data in bytes. Set to 0 if data storage is not used.
};

// Possible return values of rebound_integrate
enum REB_STATUS {
    // Any status less than SINGLE_STEP get incremented once every timestep until SINGLE_STEP is reached.
    REB_STATUS_SINGLE_STEP = -10, // Performing a single step, then switching to PAUSED.
    REB_STATUS_SCREENSHOT_READY=-5,// Screenshot is ready, send back, then finish integration
    REB_STATUS_SCREENSHOT = -4,   // Pause until visualization has taken a screenshot.
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

// Main REBOUND Simulation structure
// Note: only variables that should be accessed by users are documented here.
struct reb_simulation {
    double  t;                      // Current simulation time. Default: 0.
    double  G;                      // Gravitational constant. Default: 1.
    double  softening;              // Gravitational softening. Default: 0.
    double  dt;                     // Timestep. Default: 0.001.
    double  dt_last_done;           // Last successful timestep.
    uint64_t steps_done;            // Number of timesteps done.

    // Main particles array
    size_t  N;                      // Number of particles (includes variational particles). Default: 0.
    size_t  N_allocated;            // Current maximum space allocated in the particles array on this node. 
    struct reb_particle* particles; // Main particle array with active, variational, and test particles.
    
    // Collision routines and some integrators can operate on only some particles.
    // This map defines which ones. Set to NULL to operate on all particles.
    size_t N_map;                   // If map is not NULL, use map to operate on only N_map particles
    size_t* map;                    // Note: memory not owned and thus not freed by simulation.

    // Variational particles array
    size_t  N_var;                  // Number of variational particles. Default 0.
    struct reb_particle* particles_var;

    size_t  N_var_config; 
    struct  reb_variational_configuration* var_config;   // Configuration structs. These contain details on variational particles. 

    int     var_rescale_warning;    
    size_t  N_active;               // Number of active (i.e. not test-particle) particles. Default: -1 (all particles are active). 
    int     testparticle_type;      // 0 (default): active particles do not feel test-particles, 1: active particles feel test-particles
    int     testparticle_hidewarnings;
    char**  name_list;              // List of names used to identify particles. Managed by REBOUND. Do not directly edit/access.
    size_t  N_name_list;            // Number of entries in reb_name_list.
    struct reb_name_hash_item*    name_hash_table;        // Internal use only. Speeds up name search.
    struct reb_vec3d* gravity_cs; 
    size_t N_allocated_gravity_cs;
    struct reb_treecell** tree_root;// Pointer to the temporary tree structure.
    double opening_angle2;          // Opening angle for tree-based gravity calculation. Default 0.25.
    enum REB_STATUS status;         // Current simulation status
    int     exact_finish_time;      // 1 (default): integrate exactly to the time requested and adjust timestep if needed, 0: may overshoot by one timestep

    int force_is_velocity_dependent; // 0 (default): force only depends on position, 1: force also depends on velocities
    enum {
        REB_GRAVITY_IGNORE_TERMS_NONE = 0,              // Include all pairwise interactions in gravity calculations.
        REB_GRAVITY_IGNORE_TERMS_BETWEEN_0_AND_1 = 1,   // Ignore interactions between particles 0 and 1 in gravity calculations (used for WHFast)
        REB_GRAVITY_IGNORE_TERMS_INVOLVING_0 = 2,       // Ignore all interactions with particle 0 in gravity calculations (used for WHFast)
    } gravity_ignore_terms;         // This flag determines which interactions are included in the gravity calculation.
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
    uint32_t python_unit_l;         // Only used when working with units in python.
    uint32_t python_unit_m;         // Only used when working with units in python.
    uint32_t python_unit_t;         // Only used when working with units in python.

    // Simulation domain and ghost boxes 
    double  root_size;              // Size of a root box. 
    size_t  N_root_x;               // Number of ghost boxes in x direction. Do not change manually.
    size_t  N_root_y;
    size_t  N_root_z;
    int  N_ghost_x;              // Number of ghost boxes in x direction.
    int  N_ghost_y;
    int  N_ghost_z;

    // MPI Parallelization
#ifdef MPI
    int    mpi_id;                              // Unique id of this node (starting at 0). Used for MPI only.
    int    mpi_num;                             // Number of MPI nodes. Used for MPI only.
    struct reb_particle** particles_send;       // Send buffer for particles. There is one buffer per node.
    int*   N_particles_send;                    // Current length of particle send buffer.
    int*   N_particles_send_max;                // Maximal length of particle send buffer before realloc() is needed.
    struct reb_particle** particles_recv;       // Receive buffer for particles. There is one buffer per node.
    int*   N_particles_recv;                    // Current length of particle receive buffer.
    int*   N_particles_recv_max;                // Maximal length of particle receive buffer before realloc() is needed. */

    struct reb_treecell** tree_essential_send;  // Send buffer for cells. There is one buffer per node.
    int*   N_tree_essential_send;               // Current length of cell send buffer.
    int*   N_tree_essential_send_max;           // Maximal length of cell send buffer before realloc() is needed.
    struct reb_treecell** tree_essential_recv;  // Receive buffer for cells. There is one buffer per node.
    int*   N_tree_essential_recv;               // Current length of cell receive buffer.
    int*   N_tree_essential_recv_max;           // Maximal length of cell receive buffer before realloc() is needed.
#endif // MPI

    struct reb_collision* collisions;       // Array of current collisions. Do not change manually
    size_t N_allocated_collisions;
    unsigned int collisions_N;              // Number of collisions found during last collision search.
    double minimum_collision_velocity;      // Ensure relative velocity during collisions is at least this much (to avoid particles sinking into each other)
    double collisions_plog;                 // Keeping track of momentum transfer in collisions (for ring simulations)
    int64_t collisions_log_n;               // Cumulative number of collisions in entire simulation.

    // MEGNO Chaos indicator. These variables should not be accessed directly. Use functions provided instead.
    int calculate_megno;    // Do not change manually. Internal flag that determines if megno is calculated (default=0, but megno_init() sets it to the index of variational particles used for megno)
    double megno_Ys;        // Running megno sum (internal use)
    double megno_Yss;       // Running megno sum (internal use)
    double megno_cov_Yt;    // covariance of MEGNO Y and t
    double megno_var_t;     // variance of t 
    double megno_mean_t;    // mean of t
    double megno_mean_Y;    // mean of MEGNO Y
    double megno_initial_t; // Time when MENGO was initialized
    int64_t megno_n;        // number of covariance updates

    unsigned int rand_seed; // seed for random number generator, used by MEGNO and other random number generators in REBOUND.

    // Simulationarchive. These variables should not be accessed directly. Use functions provided instead. 
    int    simulationarchive_version;               // Version of the SA binary format (1,2=legacy format, 3=modern format)
    double simulationarchive_auto_interval;         // Current sampling cadence, in code units
    double simulationarchive_auto_walltime;         // Current sampling cadence, in wall time
    uint64_t simulationarchive_auto_step;           // Current sampling cadence, in time steps
    double simulationarchive_next;                  // Next output time (simulation time or wall time, depending on whether auto_interval or auto_walltime is set)
    uint64_t simulationarchive_next_step;           // Next output step (only used if auto_steps is set)
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
        REB_INTEGRATOR_TRACE = 25,      // TRACE integrator (Lu, Hernandez and Rein 2024)
        REB_INTEGRATOR_CUSTOM = 26,     // Custom, user-provided integrator.
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
        REB_GRAVITY_JACOBI = 5,         // Special gravity routine which includes the Jacobi terms for WH integrators 
        REB_GRAVITY_CUSTOM = 7,         // Custom, user or integrator provided gravity routine
    } gravity;

    void (*gravity_custom) (struct reb_simulation* const r);  // Used with REB_GRAVITY_CUSTOM

    // Datastructures for integrators
    struct reb_integrator ri_custom;                // Function pointers and data for REB_INTEGRATOR_CUSTOM
    struct reb_integrator_sei ri_sei;               // The SEI struct 
    struct reb_integrator_leapfrog ri_leapfrog;     // The Leapfrog struct 
    struct reb_integrator_whfast ri_whfast;         // The WHFast struct 
    struct reb_integrator_whfast512 ri_whfast512;   // The WHFast512 struct 
    struct reb_integrator_saba ri_saba;             // The SABA struct 
    struct reb_integrator_ias15 ri_ias15;           // The IAS15 struct
    struct reb_integrator_mercurius ri_mercurius;   // The MERCURIUS struct
    struct reb_integrator_trace ri_trace;           // The TRACE struct
    struct reb_integrator_janus ri_janus;           // The JANUS struct 
    struct reb_integrator_eos ri_eos;               // The EOS struct 
    struct reb_integrator_bs ri_bs;                 // The BS struct
    
    // Internal callback functions  
    void (*did_add_particle)(struct reb_simulation* r); // This callback function gets called after a particle was added. Used by some integrators to handle particle additions during timesteps.
    void (*will_remove_particle)(struct reb_simulation* r, size_t pt); // This callback function gets called before a particle will be removed. Used by some integrators to handle particle removals during timesteps.

    // ODEs. Do not access these variables directly. Use functions provided instead.
    struct reb_ode** odes;      // all ode sets (includes nbody if BS is set as integrator)
    size_t N_odes;                 // number of ode sets
    size_t N_allocated_odes;
    int ode_warnings;

    // Callback functions
    void (*additional_forces) (struct reb_simulation* const r);             // Implement any additional (non-gravitational) forces here.
    void (*pre_timestep_modifications) (struct reb_simulation* const r);    // Executed just before eaach timestep. Used by REBOUNDx.
    void (*post_timestep_modifications) (struct reb_simulation* const r);   // Executed just after each timestep. Used by REBOUNDx.
    void (*heartbeat) (struct reb_simulation* r);                           // Executed at each timestep once. Use this to do extra output/work during a simulation.
    int (*key_callback) (struct reb_simulation* r, int key);                // Used when SERVER or OPENGL visualization is on. Gets called after completed timestep and if a key has been pressed. Return 1 if you want to skip default key commands.
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v);  // Allows for a velocity dependent coefficient of restitution (used for ring simulations)
    enum REB_COLLISION_RESOLVE_OUTCOME (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);        // Determines what happens when two particles collide.
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
// Reset all integrator variables. Use this to clear any precalculated data, e.g. after changing the timestep.
DLLEXPORT void reb_simulation_integrators_reset(struct reb_simulation* r);
// Make a deep copy of simulation.
DLLEXPORT struct reb_simulation* reb_simulation_copy(struct reb_simulation* r);
// Compare r1 to r2. If exactly equal then 0 is returned, otherwise 1. If output_option=1, then difference is also printed on screen.
DLLEXPORT int reb_simulation_diff(struct reb_simulation* r1, struct reb_simulation* r2, int output_option);

// Start webserver for visualization. Returns 0 on success.
DLLEXPORT int reb_simulation_start_server(struct reb_simulation* r, int port);
// Stop webserver.
DLLEXPORT void reb_simulation_stop_server(struct reb_simulation* r);


// Errors, warnings

// For fatal errors only. Print out an error message, then exit immediately and kill the process. Does no clean up memory.
DLLEXPORT void reb_exit(const char* const msg);
// Stop current integration in a nice way. Can be called from within heartbeat function.
DLLEXPORT void reb_simulation_stop(struct reb_simulation* const r);
// Print or store an info message, then continue.
DLLEXPORT void reb_simulation_info(struct reb_simulation* const r, const char* const msg);
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
// Write a screenshot of the current simulation to a file. Requires that a server was started with reb_simulation_start_server() and one client web browser is connected. Returns 0 if successful.
DLLEXPORT int reb_simulation_output_screenshot(struct reb_simulation* r, const char* filename);


// Timestepping

// Advance simulation by N_steps timesteps.
DLLEXPORT void reb_simulation_steps(struct reb_simulation* const r, size_t N_steps);
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
// Same as above but subtract r2 from r component wise.
DLLEXPORT int reb_simulation_isub(struct reb_simulation* r, struct reb_simulation* r2);
// Updates the acceleration of all particles but does not perform a step. Used  by REBOUNDx.
DLLEXPORT void reb_simulation_update_acceleration(struct reb_simulation* r);


// Diangnostic functions

// Return the sum of potential and kinetic energy
DLLEXPORT double reb_simulation_energy(struct reb_simulation* const r);
// Returns the angular momentum.
DLLEXPORT struct reb_vec3d reb_simulation_angular_momentum(const struct reb_simulation* const r);
// Returns the center of mass of a simulation.
DLLEXPORT struct reb_particle reb_simulation_com(struct reb_simulation* r);
// Returns the center of mass of two particles.
DLLEXPORT struct reb_particle reb_particle_com_of_pair(struct reb_particle p1, struct reb_particle p2);
// Returns the center of mass of particles in the simulation within a given range.
DLLEXPORT struct reb_particle reb_simulation_com_range(struct reb_simulation* r, size_t first, size_t last);
// Returns the gravitational timescale as calculated in Pham, Rein, Spiegel (2023). Useful for setting the initial IAS15 timestep.
DLLEXPORT double reb_integrator_ias15_timescale(struct reb_simulation* r);

// Functions to add and initialize particles

// Use this function to add particle pt to simulation r.
DLLEXPORT void reb_simulation_add(struct reb_simulation* const r, struct reb_particle pt);
// Use this function to add a particle to a simulation using orbital parameters or cartesian coordinates. See examples on usage.
DLLEXPORT void reb_simulation_add_fmt(struct reb_simulation* r, const char* fmt, ...);
// This function sets up a Plummer sphere, N=number of particles, M=total mass, R=characteristic radius. Particles get added to simulation.
DLLEXPORT void reb_simulation_add_plummer(struct reb_simulation* r, size_t _N, double M, double R); 
// Returns a particle given a set of orbital parameters. Also sets err to error code if initialization failed.
DLLEXPORT struct reb_particle reb_particle_from_orbit_err(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f, int* err);
// Same as above but without error code.
DLLEXPORT struct reb_particle reb_particle_from_orbit(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f);
// Returns a particle given a set of Pal orbital parameters.
DLLEXPORT struct reb_particle reb_particle_from_pal(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy);
// Returns a reb_particle structure with fields/name/ptrs initialized to nan/0/NULL. 
DLLEXPORT struct reb_particle reb_particle_nan(void);


// Functions to access, remove, and operate on particles

// Remove all particles
DLLEXPORT void reb_simulation_remove_all_particles(struct reb_simulation* const r);
// Remove one particle. Remaining particle will keep order. Returns 0 if successful.
DLLEXPORT int reb_simulation_remove_particle(struct reb_simulation* const r, size_t index);
// Remove one particle by name. If multiple particles share the name, one particle will be remove. Which one is undetermined. Returns 0 if successful.
DLLEXPORT int reb_simulation_remove_particle_by_name(struct reb_simulation* r, const char* const name); 
// Returns a particle's index in the simulation given a pointer to the particle. Returns -1 if not found. 
DLLEXPORT int reb_simulation_particle_index(struct reb_particle* p); 
// Returns a variational particle's index in the simulation given a pointer to the particle. Returns -1 if not found. 
DLLEXPORT int reb_simulation_particle_var_index(struct reb_particle* p); 
// Subtract particle p2 from p1
DLLEXPORT void reb_particle_isub(struct reb_particle* p1, struct reb_particle* p2);
// Add particle p2 to p1
DLLEXPORT void reb_particle_iadd(struct reb_particle* p1, struct reb_particle* p2);
// Multiply x,y,z,vx,vy,vz,m of p1 with given value.
DLLEXPORT void reb_particle_imul(struct reb_particle* p1, double value);
// Return the distance between particle p1 and p2
DLLEXPORT double reb_particle_distance(struct reb_particle* p1, struct reb_particle* p2);
// Compares two particles, ignoring pointers. Returns 1 if particles differ, 0 if they are exactly equal.
DLLEXPORT int reb_particle_cmp(struct reb_particle p1, struct reb_particle p2); 
// Advances one particle forward in a Keplerian orbit for time dt. mu is the gravitational parameter, G*(m+M). Set r=NULL unless variational particles are used. Returns 0 on success, 1 if timestep is too large. 
DLLEXPORT int reb_integrator_whfast_kepler_solver(struct reb_particle* const restrict p, double mu, double dt, const struct reb_simulation* const r);
// Sets a particle's name. This function should be used instead of directly setting the name in the particle's structure as it
// registers the name, allowing for faster lookup and storing of name in binary files.
DLLEXPORT void reb_particle_set_name(struct reb_particle* p, const char* const name);
/// Register a string. This copies the string, returns the new pointer. REBOUND will manage the memory of the copy. Used by REBOUNDx.
DLLEXPORT const char* reb_simulation_register_name(struct reb_simulation* r, const char* const name);
// Returns a pointer to a particle given its name. Returns NULL if not found.
DLLEXPORT struct reb_particle* reb_simulation_get_particle_by_name(struct reb_simulation* r, const char* const name);
#ifdef MPI
// When MPI is used, particles cannot be accessed by name. Need to use id instead.
DLLEXPORT struct reb_particle reb_simulation_particle_by_id(struct reb_simulation* const r, size_t id);
#endif // MPI



// Chaos indicators

// Turn on MEGNO/Lyapunov calculation. Uses random seen in simulation.
DLLEXPORT void reb_simulation_init_megno(struct reb_simulation* const r);
// Same as above but used given random seed. Useful to reproduce same results every time.
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

// Built in trace switching functions

DLLEXPORT int reb_integrator_trace_switch_peri_default(struct reb_simulation* const r, const size_t j);
DLLEXPORT int reb_integrator_trace_switch_peri_none(struct reb_simulation* const r, const size_t j);
DLLEXPORT int reb_integrator_trace_switch_default(struct reb_simulation* const r, const size_t i, const size_t j);

// Built in collision resolve functions

DLLEXPORT enum REB_COLLISION_RESOLVE_OUTCOME reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c); // halts a simulation when a collision occurs
DLLEXPORT enum REB_COLLISION_RESOLVE_OUTCOME reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);
DLLEXPORT enum REB_COLLISION_RESOLVE_OUTCOME reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);

// Random sampling - These functions only use the simulation object for a seed. If r=NULL, time and PID are used as a seed.

DLLEXPORT double reb_random_uniform(struct reb_simulation* r, double min, double max);
DLLEXPORT double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
DLLEXPORT double reb_random_normal(struct reb_simulation* r, double variance);
DLLEXPORT double reb_random_rayleigh(struct reb_simulation* r, double sigma);


// Miscellaneous functions

// Returns the angle f wrapped in the interval from 0 to 2*pi
DLLEXPORT double reb_mod2pi(double f);
// True anomaly for a given eccentricity and mean anomaly
DLLEXPORT double reb_M_to_f(double e, double M);
// True anomaly for a given eccentricity and eccentric anomaly
DLLEXPORT double reb_E_to_f(double e, double M);
// Eccentric anomaly for a given eccentricity and mean anomaly
DLLEXPORT double reb_M_to_E(double e, double M);
// Returns the hash for a given string. 
DLLEXPORT uint32_t reb_hash(const char* c);


// Simulationarchive

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
// Same as above but returns sets an error code if not successful.
DLLEXPORT struct reb_orbit reb_orbit_from_particle_err(double G, struct reb_particle p, struct reb_particle primary, int* err);


// ODE functions

// Defines one Ordinary Differential Equation (ODE) so that it can be integrated with REBOUND
struct reb_ode{
    unsigned int length;        // number of components / dimension
    double* y;                  // Pointer to current state 
    unsigned int needs_nbody;   // 1: ODE needs N-body particles to calculate RHS
    void* ref;                  // Optional pointer to any additional data needed for derivative calculation
    void (*derivatives)(struct reb_ode* const ode, double* const yDot, const double* const y, const double t);  // Function pointer to right hand side of ODE
    void (*getscale)(struct reb_ode* const ode, const double* const y0, const double* const y1);                // Function pointer, sets scales for components (optional) 
    void (*pre_timestep)(struct reb_ode* const ode, const double* const y0);                                    // Function pointer, gets called just before the ODE integration (optional)
    void (*post_timestep)(struct reb_ode* const ode, const double* const y0);                                   // Function pointer, gets called just after the ODE integration (optional)

    // Internal use
    size_t N_allocated;   
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
DLLEXPORT int reb_simulation_add_variation_2nd_order(struct reb_simulation* const r, int testparticle, size_t index_1st_order_a, size_t index_1st_order_b);

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
//   implemented as separate functions. 
//   The following variables are supported: a, e, inc, f, omega, Omega, h, k, ix, iy and m (mass). 
// The functions return the derivative as a particle structure. Each structure element is a derivative.
// The parameter po is the original particle for which the derivative is to be calculated.
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

// Functions for Frequency Analysis, MFT, FMFT
enum REB_FREQUENCY_ANALYSIS_TYPE {
    REB_FREQUENCY_ANALYSIS_MFT = 0,
    REB_FREQUENCY_ANALYSIS_FMFT = 1,
    REB_FREQUENCY_ANALYSIS_FMFT2 = 2,
};
// Returns 0 on success
DLLEXPORT int reb_frequency_analysis(double *output, size_t nfreq, double minfreq, double maxfreq, enum REB_FREQUENCY_ANALYSIS_TYPE type, double *input, size_t ndata);

// Functions to convert between coordinate systems

// Rotations

// Rotation struct (implemented as a quaternion)
struct reb_rotation {
    double ix;
    double iy;
    double iz;
    double r;
};

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


// The following structures are related to OpenGL/WebGL visualization.

// Generic 3d vector
struct reb_vec3df {
    float x,y,z;
};

// Generic 4d matrix (single precision)
struct reb_mat4df {
    float m[16];
};

// Structures and functions used for programmatic animations.
// Current display orientation, settings, etc
struct reb_display_settings {
    struct reb_mat4df view;
    int spheres;                    // Switches between point sprite and real spheres.
    int pause;                      // Pauses visualization, but keep simulation running
    int wire;                       // Shows/hides orbit wires.
    unsigned int breadcrumbs;       // Number of past particle positions.
    int onscreentext;               // Shows/hides onscreen text.
    int onscreenhelp;               // Shows/hides onscreen help.
    int multisample;                // Turn off/on multisampling.
    int ghostboxes;                 // Shows/hides ghost boxes.
    int reference;                  // reb_particle used as a reference for centering.
};

// Display settings initialization. Overwrites user interactions.
DLLEXPORT void reb_simulation_add_display_settings(struct reb_simulation* r);

// Matrix/vector methods for single precision operations. Used for graphics only
DLLEXPORT struct reb_mat4df reb_mat4df_identity();
DLLEXPORT struct reb_mat4df reb_mat4df_scale(struct reb_mat4df m, float x, float y, float z);
DLLEXPORT void reb_mat4df_print(struct reb_mat4df m);
DLLEXPORT int reb_mat4df_eq(struct reb_mat4df A, struct reb_mat4df B);
DLLEXPORT struct reb_vec3df reb_mat4df_get_scale(struct reb_mat4df m);
DLLEXPORT struct reb_mat4df reb_mat4df_translate(struct reb_mat4df m, float x, float y, float z);
DLLEXPORT struct reb_mat4df reb_mat4df_multiply(struct reb_mat4df A, struct reb_mat4df B);
DLLEXPORT struct reb_mat4df reb_rotation_to_mat4df(struct reb_rotation A);
DLLEXPORT struct reb_mat4df reb_mat4df_ortho(float l, float r, float b, float t, float n, float f);


#endif // _MAIN_H
