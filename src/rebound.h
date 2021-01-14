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

#ifdef __cplusplus

// At least GCC and clang support the restrict keyword as an extension.
#if defined(__GNUC__) || defined(__clang__)

#define REBOUND_RESTRICT __restrict__

#else

// For other compilers, we disable it.
#define REBOUND_RESTRICT

#endif

extern "C" {

#else

#define REBOUND_RESTRICT restrict

#endif

#include <inttypes.h>
#include <stdint.h>
#include <sys/time.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#ifndef M_PI
// Make sure M_PI is defined. 
#define M_PI           3.14159265358979323846       ///< The mathematical constant pi.
#endif
#ifdef MPI
#include "mpi.h"
#endif // MPI
#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif // GITHASH

extern const char* reb_build_str;   ///< Date and time build string.
extern const char* reb_version_str; ///< Version string.
extern const char* reb_githash_str; ///< Current git hash.
extern const char* reb_logo[26];    ///< Logo of rebound. 
extern volatile sig_atomic_t reb_sigint;  ///< Graceful global interrupt handler 

// Forward declarations
struct reb_simulation;
struct reb_display_data;
struct reb_treecell;

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
    uint32_t hash;      ///< hash to identify particle.
    void* ap;           ///< Functionality for externally adding additional properties to particles.
    struct reb_simulation* sim; ///< Pointer to the parent simulation.
};

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
    double* REBOUND_RESTRICT p0; ///< 0 substep
    double* REBOUND_RESTRICT p1; ///< 1 substep
    double* REBOUND_RESTRICT p2; ///< 2 substep
    double* REBOUND_RESTRICT p3; ///< 3 substep
    double* REBOUND_RESTRICT p4; ///< 4 substep
    double* REBOUND_RESTRICT p5; ///< 5 substep
    double* REBOUND_RESTRICT p6; ///< 6 substep
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
 * @defgroup IntegratorStructs Configuration structures for integrators
 * @{
*/
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
     * @brief This parameter controls the order of operations.
     * @details The order of floating point operations can affect the long term energy error. In the original version of IAS15, 
     * the order of floating point operations leads to a linear energy error growth when a fixed timestep is used (Hernandez + Holman 2019).
     * By default, this parameter is now set to 1 which uses a different order and thus avoids this issue. Note that using IAS15 as a 
     * non-adaptive integrator is rarely beneficial. This flag only exists to allow for backwards compatibility. All new simulations should 
     * use neworder=1.
     **/
    unsigned int neworder;

    
    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    /**
     * @brief Counter how many times the iteration did not converge. 
     */
    unsigned long iterations_max_exceeded;



    int allocatedN;             ///< Size of allocated arrays.

    double* REBOUND_RESTRICT at;            ///< Temporary buffer for acceleration
    double* REBOUND_RESTRICT x0;            ///<                      position (used for initial values at h=0)
    double* REBOUND_RESTRICT v0;            ///<                      velocity
    double* REBOUND_RESTRICT a0;            ///<                      acceleration
    double* REBOUND_RESTRICT csx;           ///<                      compensated summation for x
    double* REBOUND_RESTRICT csv;           ///<                      compensated summation for v
    double* REBOUND_RESTRICT csa0;          ///<                      compensated summation for a

    struct reb_dp7 g;
    struct reb_dp7 b;
    struct reb_dp7 csb;         ///< Compensated summation for b
    struct reb_dp7 e;

    // The following values are used for resetting the b and e coefficients if a timestep gets rejected
    struct reb_dp7 br;
    struct reb_dp7 er;

    int* map;               // map to particles (identity map for non-mercurius simulations)
    int map_allocated_N;    // allocated size for map
    /**
     * @endcond
     */

};

/**
 * @brief This structure contains variables and pointer used by the MERCURIUS integrator.
 */
struct reb_simulation_integrator_mercurius {
   /**
    * @brief This is a function pointer to the force switching function used.
    * @details If NULL (the default), the MERCURY switching function will be used.
    * The argument d is the distance between two particles.
    * The argument dcrit is the maximum of the critical distances of the two particles.
    * The return value is a scalar between 0 and 1. If it always returns 1, then the
    * integrator becomes the standard Wisdom-Holman integrator.
    */
    double (*L) (const struct reb_simulation* const r, double d, double dcrit);  
    
    /** 
     * @brief Switching distance in units of the Hill radius 
     * @brief The switching distances for particles are calculated automatically
     * based on multiple criteria. One criterion calculates the Hill radius of 
     * particles and then multiplies it with the hillfac variable. 
     */ 
    double hillfac;        
    
    /** 
     * @brief Setting this flag to one will recalculate heliocentric coordinates from the particle structure at the beginning of the next timestep. 
     * @details After one timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_coordinates_this_timestep;

    /** 
     * @brief Setting this flag to one will recalculate the critical switchover 
     * distances dcrit at the beginning of the next timestep. 
     * @details After one timestep, the flag gets set back to 0. 
     * If you want to recalculate dcrit at every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_dcrit_this_timestep;

    /**
     * @brief If this flag is set (the default), the integrator will 
     * recalculate heliocentric coordinates and synchronize after
     * every timestep, to avoid problems with outputs or particle modifications
     * between timesteps. 
     * @details Setting it to 0 will result in a speedup, but care
     * must be taken to synchronize and recalculate coordinates when needed.
     */
    unsigned int safe_mode;
    
    unsigned int is_synchronized;   ///< Flag to determine if current particle structure is synchronized
    unsigned int mode;              ///< Internal. 0 if WH is operating, 1 if IAS15 is operating.
    unsigned int encounterN;        ///< Number of particles currently having an encounter
    unsigned int encounterNactive;  ///< Number of particles currently having an encounter
    unsigned int allocatedN;        ///< Current size of allocated internal arrays
    unsigned int allocatedN_additionalforces;        ///< Current size of allocated internal particles_backup_additionalforces array
    unsigned int dcrit_allocatedN;  ///< Current size of dcrit arrays
    double* dcrit;                  ///< Switching radii for particles
    struct reb_particle* REBOUND_RESTRICT particles_backup;     ///< Internal array, contains coordinates before Kepler step for encounter prediction
    struct reb_particle* REBOUND_RESTRICT particles_backup_additionalforces;     ///< Internal array, contains coordinates before Kepler step for encounter prediction
    int* encounter_map;             ///< Map to represent which particles are integrated with ias15
    struct reb_vec3d com_pos;       ///< Used internally to keep track of the centre of mass during the timestep
    struct reb_vec3d com_vel;       ///< Used internally to keep track of the centre of mass during the timestep
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
 * @brief This structure contains variables used by the SABA integrator.
 */
struct reb_simulation_integrator_saba {
    /**
     * @brief SABA type.
     * @details Available types include: SABA1, SABA2, SABA3, SABA4, SABACM1, 
     * SABACM2, SABACM3, SABACM4, SABACL1, SABACL2, SABACL3, SABACL4, 
     * SABA(10,4), SABA(8,6,4), SABA(10,6,4), SABAH(8,4,4), SABAH(8,6,4), 
     * and SABAH(10,6,4).
     */
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
    } type;
    unsigned int safe_mode;       ///< Safe_mode has the same functionality as in WHFast.
    unsigned int is_synchronized; ///< Flag to determine if current particle structure is synchronized
    /**
     * @brief Flaf that determines if the inertial coordinates generated are discarded in subsequent timesteps (Jacobi coordinates are used instead).
     * @details Danger zone! Only use this flag if you are absolutely sure
     * what you are doing. This is intended for
     * simulation which have to be reproducible on a bit by bit basis.
     */
    unsigned int keep_unsynchronized;
};

/**
 * @brief This structure contains variables used by the WHFast integrator.
 */
struct reb_simulation_integrator_whfast {
    /**
     * @brief This variable turns on/off different first symplectic correctors for WHFast.
     * @details These correctors remove terms of order O(eps*dt^2) 
     * - 0 (default): turns off all first correctors
     * - 3: uses third order (two-stage) first corrector 
     * - 5: uses fifth order (four-stage) first corrector 
     * - 7: uses seventh order (six-stage) first corrector 
     * - 11: uses eleventh order (ten-stage) first corrector 
     * - 17: uses 17th order (16-stage) first corrector 
     */
    unsigned int corrector;
    
    /**
     * @brief This variable turns on/off the second symplectic correctors for WHFast.
     * @details 
     * - 0 (default): turns off second correctors
     * - 1: uses second corrector 
     */
    unsigned int corrector2;
    
    /**
     * @brief This variable determines the kernel of the WHFast integrator.
     * @details 
     * - 0 (default): Uses a standard WH kick step 
     * - 1: uses the exact modified kick (for Newtonian gravity) 
     * - 2: uses the composition kernel  
     * - 3: uses the lazy implementer's modified kick   
     */
    enum {
        REB_WHFAST_KERNEL_DEFAULT = 0,
        REB_WHFAST_KERNEL_MODIFIEDKICK = 1,
        REB_WHFAST_KERNEL_COMPOSITION = 2,
        REB_WHFAST_KERNEL_LAZY = 3,
    }kernel;
    
    
    /**
     * @brief Chooses the coordinate system for the WHFast algorithm. Default is Jacobi Coordinates.
     */
    enum {
        REB_WHFAST_COORDINATES_JACOBI = 0,                      ///< Jacobi coordinates (default)
        REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC = 1,      ///< Democratic Heliocentric coordinates
        REB_WHFAST_COORDINATES_WHDS = 2,                        ///< WHDS coordinates (Hernandez and Dehnen, 2017)
        } coordinates;

    /** 
     * @brief Setting this flag to one will recalculate Jacobi/heliocentric coordinates from the particle structure in the next timestep. 
     * @details After the timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_coordinates_this_timestep;

    /**
     * @brief If this flag is set (the default), whfast will recalculate 
     * jacobi/heliocentric coordinates and synchronize
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
     * @brief Jacobi/heliocentric coordinates
     * @details This array contains the Jacobi/heliocentric
     * coordinates of all particles.
     * It is automatically filled and updated by WHfast.
     * Access this array with caution.
     */
    struct reb_particle* REBOUND_RESTRICT p_jh;
    
    /**
     * @brief Internal temporary array used for lazy implementer's kernel method
     */
    struct reb_particle* REBOUND_RESTRICT p_temp;
    
    /**
     * @brief Generate inertial coordinates at the end of the integration, but do not change the Jacobi/heliocentric coordinates
     * @details Danger zone! Only use this flag if you are absolutely sure
     * what you are doing. This is intended for
     * simulation which have to be reproducible on a bit by bit basis.
     */
    unsigned int keep_unsynchronized;

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */

    unsigned int is_synchronized;   ///< Flag to determine if current particle structure is synchronized
    unsigned int allocated_N;       ///< Space allocated in p_jh array
    unsigned int allocated_Ntemp;   ///< Space allocated in p_temp array
    unsigned int timestep_warning;  ///< Counter of timestep warnings
    unsigned int recalculate_coordinates_but_not_synchronized_warning;   ///< Counter of Jacobi synchronization errors
    /**
     * @endcond
     */
};

/**
 * @brief Available operator splitting methods for phi0 and phi1 in EOS integrators.
 */
enum REB_EOS_TYPE {
    REB_EOS_LF = 0x00,      // 2nd order, standard leap-frog
    REB_EOS_LF4 = 0x01,     // 4th order, three function evaluations
    REB_EOS_LF6 = 0x02,     // 6th order, nine function evaluations
    REB_EOS_LF8 = 0x03,     // 8th order, seventeen funtion evaluations, see Blanes & Casa (2016), p91
    REB_EOS_LF4_2 = 0x04,   // generalized order (4,2), two force evaluations, McLachlan 1995
    REB_EOS_LF8_6_4= 0x05,  // generalized order (8,6,4), seven force evaluations
    REB_EOS_PLF7_6_4= 0x06, // generalized order (7,6,4), three force evaluations, pre- and post-processors
    REB_EOS_PMLF4 = 0x07,   // 4th order, one modified force evaluation, pre- and post-processors, Blanes et al. (1999)
    REB_EOS_PMLF6 = 0x08,   // 6th order, three modified force evaluations, pre- and post-processors, Blanes et al. (1999)
};

/**
 * @brief This structure contains variables used by the EOS integrator family.
 */
struct reb_simulation_integrator_eos {
    enum REB_EOS_TYPE phi0;         ///< Outer operator splitting scheme
    enum REB_EOS_TYPE phi1;         ///< Inner operator splitting scheme
    unsigned int n;                 ///

    unsigned int safe_mode;         ///< If set to 0, always combine drift steps at the beginning and end of phi0. If set to 1, n needs to be bigger than 1.
    unsigned int is_synchronized;   ///< Flag to indicate if the drift step at the end of the last timestep has been taken.
};



/**
 * @cond PRIVATE
 * Internal data structures below. Nothing to be changed by the user.
 */
#define REB_PARTICLE_INT_TYPE int64_t
/**
 * @brief Integer positions and velocities for particles. Used in JANUS integrator. 
 */
struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};
/**
 * @endcond
 */

struct reb_simulation_integrator_janus {
    /**
     * @brief Scale of the problem. Positions get divided by this number before the conversion to an integer. 
     */
    double scale_pos;
    /**
     * @brief Scale of the problem. Velocities get divided by this number before the conversion to an integer. 
     */
    double scale_vel;
    /**
     * @brief Order of the scheme. Default is 6. 
     */
    unsigned int order; //TODO needs input/output
    /**
     * @brief If this flag is set, then janus will recalculate integer coordinates at
     * the next timestep.
     */
    unsigned int recalculate_integer_coordinates_this_timestep;
    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    struct reb_particle_int* REBOUND_RESTRICT p_int;    ///< Integer particle pos/vel
    unsigned int allocated_N;                   ///< Space allocated in arrays
    /**
     * @endcond
     */
};

/**
 * @defgroup MiscRebStructs Miscellaneous REBOUND structures
 * @{
*/

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
 * @brief Enumeration describing the return status of rebound_integrate
 */
enum REB_STATUS {
    REB_RUNNING_PAUSED = -3,    ///< Simulation is paused by visualization.
    REB_RUNNING_LAST_STEP = -2, ///< Current timestep is the last one. Needed to ensure that t=tmax exactly.
    REB_RUNNING = -1,           ///< Simulation is current running, no error occurred.
    REB_EXIT_SUCCESS = 0,       ///< Integration finished successfully.
    REB_EXIT_ERROR = 1,         ///< A generic error occurred and the integration was not successful.
    REB_EXIT_NOPARTICLES = 2,   ///< The integration ends early because no particles are left in the simulation.
    REB_EXIT_ENCOUNTER = 3,     ///< The integration ends early because two particles had a close encounter (see exit_min_distance)
    REB_EXIT_ESCAPE = 4,        ///< The integration ends early because a particle escaped (see exit_max_distance)  
    REB_EXIT_USER = 5,          ///< User caused exit, simulation did not finish successfully.
    REB_EXIT_SIGINT = 6,        ///< SIGINT received. Simulation stopped.
    REB_EXIT_COLLISION = 7,     ///< The integration ends early because two particles collided. 
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
    int index_1st_order_a;      ///< Used for 2nd order variational particles only: Index of the first order variational particle in the particles array.
    int index_1st_order_b;      ///< Used for 2nd order variational particles only: Index of the first order variational particle in the particles array.
};

/**
 * @cond PRIVATE
 * Internal data structures below. Nothing to be changed by the user.
 */

/**
 * @brief Enumeration describing the contents of a binary field. Used to read and write binary files.
 */
enum REB_BINARY_FIELD_TYPE {
    REB_BINARY_FIELD_TYPE_T = 0,
    REB_BINARY_FIELD_TYPE_G = 1,
    REB_BINARY_FIELD_TYPE_SOFTENING = 2,
    REB_BINARY_FIELD_TYPE_DT = 3,
    REB_BINARY_FIELD_TYPE_N = 4,
    REB_BINARY_FIELD_TYPE_NVAR = 5,
    REB_BINARY_FIELD_TYPE_VARCONFIGN = 6,
    REB_BINARY_FIELD_TYPE_NACTIVE = 7,
    REB_BINARY_FIELD_TYPE_TESTPARTICLETYPE = 8,
    REB_BINARY_FIELD_TYPE_HASHCTR = 9, 
    REB_BINARY_FIELD_TYPE_OPENINGANGLE2 = 10,
    REB_BINARY_FIELD_TYPE_STATUS = 11,
    REB_BINARY_FIELD_TYPE_EXACTFINISHTIME = 12,
    REB_BINARY_FIELD_TYPE_FORCEISVELOCITYDEP = 13,
    REB_BINARY_FIELD_TYPE_GRAVITYIGNORETERMS = 14,
    REB_BINARY_FIELD_TYPE_OUTPUTTIMINGLAST = 15,
    REB_BINARY_FIELD_TYPE_SAVEMESSAGES = 16,
    REB_BINARY_FIELD_TYPE_EXITMAXDISTANCE = 17,
    REB_BINARY_FIELD_TYPE_EXITMINDISTANCE = 18,
    REB_BINARY_FIELD_TYPE_USLEEP = 19,
    REB_BINARY_FIELD_TYPE_TRACKENERGYOFFSET = 20,
    REB_BINARY_FIELD_TYPE_ENERGYOFFSET = 21,
    REB_BINARY_FIELD_TYPE_BOXSIZE = 22, 
    REB_BINARY_FIELD_TYPE_BOXSIZEMAX = 23, 
    REB_BINARY_FIELD_TYPE_ROOTSIZE = 24,
    REB_BINARY_FIELD_TYPE_ROOTN = 25,
    REB_BINARY_FIELD_TYPE_ROOTNX = 26, 
    REB_BINARY_FIELD_TYPE_ROOTNY = 27,
    REB_BINARY_FIELD_TYPE_ROOTNZ = 28,
    REB_BINARY_FIELD_TYPE_NGHOSTX = 29,
    REB_BINARY_FIELD_TYPE_NGHOSTY = 30,
    REB_BINARY_FIELD_TYPE_NGHOSTZ = 31,
    REB_BINARY_FIELD_TYPE_COLLISIONRESOLVEKEEPSORTED = 32,
    REB_BINARY_FIELD_TYPE_MINIMUMCOLLISIONVELOCITY = 33,
    REB_BINARY_FIELD_TYPE_COLLISIONSPLOG = 34, 
    REB_BINARY_FIELD_TYPE_MAXRADIUS = 35, 
    REB_BINARY_FIELD_TYPE_COLLISIONSNLOG = 36, 
    REB_BINARY_FIELD_TYPE_CALCULATEMEGNO = 37, 
    REB_BINARY_FIELD_TYPE_MEGNOYS = 38, 
    REB_BINARY_FIELD_TYPE_MEGNOYSS = 39, 
    REB_BINARY_FIELD_TYPE_MEGNOCOVYT = 40,
    REB_BINARY_FIELD_TYPE_MEGNOVART = 41, 
    REB_BINARY_FIELD_TYPE_MEGNOMEANT = 42, 
    REB_BINARY_FIELD_TYPE_MEGNOMEANY = 43, 
    REB_BINARY_FIELD_TYPE_MEGNON = 44,
    REB_BINARY_FIELD_TYPE_SASIZEFIRST = 45,
    REB_BINARY_FIELD_TYPE_SASIZESNAPSHOT = 46,
    REB_BINARY_FIELD_TYPE_SAAUTOINTERVAL = 47,
    REB_BINARY_FIELD_TYPE_SAAUTOWALLTIME = 102,
    REB_BINARY_FIELD_TYPE_SANEXT = 48,
    REB_BINARY_FIELD_TYPE_COLLISION = 50,
    REB_BINARY_FIELD_TYPE_INTEGRATOR = 51,
    REB_BINARY_FIELD_TYPE_BOUNDARY = 52,
    REB_BINARY_FIELD_TYPE_GRAVITY = 53,
    REB_BINARY_FIELD_TYPE_SEI_OMEGA = 54,
    REB_BINARY_FIELD_TYPE_SEI_OMEGAZ = 55,
    REB_BINARY_FIELD_TYPE_SEI_LASTDT = 56,
    REB_BINARY_FIELD_TYPE_SEI_SINDT = 57,
    REB_BINARY_FIELD_TYPE_SEI_TANDT = 58,
    REB_BINARY_FIELD_TYPE_SEI_SINDTZ = 59,
    REB_BINARY_FIELD_TYPE_SEI_TANDTZ = 60,
    REB_BINARY_FIELD_TYPE_WHFAST_CORRECTOR = 61,
    REB_BINARY_FIELD_TYPE_WHFAST_RECALCJAC = 62, 
    REB_BINARY_FIELD_TYPE_WHFAST_SAFEMODE = 63,
    REB_BINARY_FIELD_TYPE_WHFAST_KEEPUNSYNC = 64,
    REB_BINARY_FIELD_TYPE_WHFAST_ISSYNCHRON = 65,
    REB_BINARY_FIELD_TYPE_WHFAST_TIMESTEPWARN = 66,
    REB_BINARY_FIELD_TYPE_IAS15_EPSILON = 69,
    REB_BINARY_FIELD_TYPE_IAS15_MINDT = 70,
    REB_BINARY_FIELD_TYPE_IAS15_EPSILONGLOBAL = 71,
    REB_BINARY_FIELD_TYPE_IAS15_ITERATIONSMAX = 72,
    REB_BINARY_FIELD_TYPE_PARTICLES = 85,
    REB_BINARY_FIELD_TYPE_VARCONFIG = 86,
    REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS = 87,
    REB_BINARY_FIELD_TYPE_IAS15_ALLOCATEDN = 88,
    REB_BINARY_FIELD_TYPE_IAS15_AT = 89,
    REB_BINARY_FIELD_TYPE_IAS15_X0 = 90,
    REB_BINARY_FIELD_TYPE_IAS15_V0 = 91,
    REB_BINARY_FIELD_TYPE_IAS15_A0 = 92,
    REB_BINARY_FIELD_TYPE_IAS15_CSX = 93,
    REB_BINARY_FIELD_TYPE_IAS15_CSV = 94,
    REB_BINARY_FIELD_TYPE_IAS15_CSA0 = 95,
    REB_BINARY_FIELD_TYPE_IAS15_G = 96,
    REB_BINARY_FIELD_TYPE_IAS15_B = 97,
    REB_BINARY_FIELD_TYPE_IAS15_CSB = 98,
    REB_BINARY_FIELD_TYPE_IAS15_E = 99,
    REB_BINARY_FIELD_TYPE_IAS15_BR = 100,
    REB_BINARY_FIELD_TYPE_IAS15_ER = 101,
    REB_BINARY_FIELD_TYPE_WHFAST_PJ = 104,
    REB_BINARY_FIELD_TYPE_VISUALIZATION = 107,
    REB_BINARY_FIELD_TYPE_JANUS_ALLOCATEDN = 110,
    REB_BINARY_FIELD_TYPE_JANUS_PINT = 112,
    REB_BINARY_FIELD_TYPE_JANUS_SCALEPOS = 113,
    REB_BINARY_FIELD_TYPE_JANUS_SCALEVEL = 114,
    REB_BINARY_FIELD_TYPE_JANUS_ORDER = 115,
    REB_BINARY_FIELD_TYPE_JANUS_RECALC = 116,
    REB_BINARY_FIELD_TYPE_WHFAST_COORDINATES = 117,
    REB_BINARY_FIELD_TYPE_MERCURIUS_HILLFAC = 118,
    REB_BINARY_FIELD_TYPE_MERCURIUS_SAFEMODE = 119,
    REB_BINARY_FIELD_TYPE_MERCURIUS_ISSYNCHRON = 120,
    REB_BINARY_FIELD_TYPE_MERCURIUS_DCRIT = 122,
    REB_BINARY_FIELD_TYPE_SAVERSION = 125,
    REB_BINARY_FIELD_TYPE_WALLTIME = 126,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_L = 130,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_M = 131,
    REB_BINARY_FIELD_TYPE_PYTHON_UNIT_T = 132,
    REB_BINARY_FIELD_TYPE_MERCURIUS_COMPOS = 133,
    REB_BINARY_FIELD_TYPE_MERCURIUS_COMVEL = 134,
    REB_BINARY_FIELD_TYPE_SAAUTOSTEP = 135,
    REB_BINARY_FIELD_TYPE_SANEXTSTEP = 136,
    REB_BINARY_FIELD_TYPE_STEPSDONE = 137,
    REB_BINARY_FIELD_TYPE_SABA_SAFEMODE = 140,
    REB_BINARY_FIELD_TYPE_SABA_ISSYNCHRON = 141,
    REB_BINARY_FIELD_TYPE_WHFAST_CORRECTOR2 = 143,
    REB_BINARY_FIELD_TYPE_WHFAST_KERNEL = 144,
    REB_BINARY_FIELD_TYPE_DTLASTDONE = 145,
    REB_BINARY_FIELD_TYPE_SABA_TYPE = 146,
    REB_BINARY_FIELD_TYPE_SABA_KEEPUNSYNC = 147,
    REB_BINARY_FIELD_TYPE_EOS_PHI0 = 148,
    REB_BINARY_FIELD_TYPE_EOS_PHI1 = 149,
    REB_BINARY_FIELD_TYPE_EOS_N = 150,
    REB_BINARY_FIELD_TYPE_EOS_SAFEMODE = 151,
    REB_BINARY_FIELD_TYPE_EOS_ISSYNCHRON = 152,
    REB_BINARY_FIELD_TYPE_IAS15_NEWORDER = 153,
    REB_BINARY_FIELD_TYPE_RAND_SEED = 154,

    REB_BINARY_FIELD_TYPE_HEADER = 1329743186,  // Corresponds to REBO (first characters of header text)
    REB_BINARY_FIELD_TYPE_SABLOB = 9998,        // SA Blob
    REB_BINARY_FIELD_TYPE_END = 9999,
};

/**
 * @brief This structure is used to save and load binary files.
 */
struct reb_binary_field {
    uint32_t type;  ///< Type of what field (enum of REB_BINARY_FIELD_TYPE)
    uint64_t size;  ///< Size in bytes of field (only counting what follows, not the binary field, itself).
};

/**
 * @brief This structure is used to save and load simulation archive files.
 */
struct reb_simulationarchive_blob {
    int32_t index;                         ///< Index of previous blob (binary file is 0, first blob is 1)
    int16_t offset_prev;                   ///< Offset to beginning of previous blob (size of previous blob).
    int16_t offset_next;                   ///< Offset to end of following blob (size of following blob).
};


/**
 * @brief This structure is used to save and load SimulationArchive files.
 * @details Everthing in this struct is handled by REBOUND itself. Users 
 * should not need to access this struct manually.
 */
struct reb_simulationarchive{
    FILE* inf;              ///< File pointer (will be kept open)
    char* filename;         ///< Filename of open file
    int version;            ///< SimulationArchive version
    long size_first;        ///< Size of first snapshot (only used for version 1)
    long size_snapshot;     ///< Size of snapshot (only used for version 1)
    double auto_interval;   ///< Interval setting used to create SA (if used)
    double auto_walltime;   ///< Walltime setting used to create SA (if used)
    unsigned long long auto_step;  ///< Steps in-between SA snapshots (if used)
    long nblobs;            ///< Total number of snapshots (including initial binary)
    uint32_t* offset;       ///< Index of offsets in file (length nblobs)
    double* t;              ///< Index of simulation times in file (length nblobs)
};

/**
 * @brief Holds a particle's hash and the particle's index in the particles array.
 * @details This structure is used for the simulation's particle_lookup_table.
 */
struct reb_hash_pointer_pair{
    uint32_t hash;
    int index;
};
/**
 * @endcond
 */
/** @} */

/**
 * @defgroup MainRebStructs Main REBOUND structures
 * @{
*/


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
    double rhill;    ///< Circular Hill radius 
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
    unsigned long long steps_done;  ///< Timesteps done
    int     N;                      ///< Current number of particles on this node. 
    int     N_var;                  ///< Total number of variational particles. Default: 0.
    int     var_config_N;           ///< Number of variational configuration structs. Default: 0.
    struct reb_variational_configuration* var_config;   ///< These configuration structs contain details on variational particles. 
    int     N_active;               ///< Number of massive particles included in force calculation (default: N). Particles with index >= N_active are considered testparticles.
    int     testparticle_type;      ///< Type of the particles with an index>=N_active. 0 means particle does not influence any other particle (default), 1 means particles with index < N_active feel testparticles (similar to MERCURY's small particles). Testparticles never feel each other.
    struct reb_hash_pointer_pair* particle_lookup_table; ///< Array of pairs that map particles' hashes to their index in the particles array.
    int     hash_ctr;               ///< Counter for number of assigned hashes to assign unique values.
    int     N_lookup;               ///< Number of entries in the particle lookup table.
    int     allocatedN_lookup;      ///< Number of lookup table entries allocated.
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
    unsigned int gravity_ignore_terms; ///< Ignore the gravity form the central object (1 for WHFast, 2 for WHFastHelio, 0 otherwise)
    double output_timing_last;      ///< Time when reb_output_timing() was called the last time. 
    unsigned long display_clock;    ///< Display clock, internal variable for timing refreshs.
    int save_messages;              ///< Set to 1 to ignore messages (used in python interface).
    char** messages;                ///< Array of strings containing last messages (only used if save_messages==1). 
    double exit_max_distance;       ///< Exit simulation if distance from origin larger than this value 
    double exit_min_distance;       ///< Exit simulation if distance from another particle smaller than this value 
    double usleep;                  ///< Wait this number of microseconds after each timestep, useful for slowing down visualization.  
    struct reb_display_data* display_data; /// < Datastructure stores visualization related data. Does not have to be modified by the user. 
    int track_energy_offset;        ///< Track energy change during collisions and ejections (default: 0).
    double energy_offset;           ///< Energy offset due to collisions and ejections (only calculated if track_energy_offset=1).
    double walltime;                ///< Walltime in seconds used by REBOUND for this simulation (integration only, not visualization, heartbeat function, etc).
    uint32_t python_unit_l;         ///< Information only used for when working with units in python.
    uint32_t python_unit_m;         ///< Information only used for when working with units in python.
    uint32_t python_unit_t;         ///< Information only used for when working with units in python.
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
    int collision_resolve_keep_sorted;      ///< Keep particles sorted if collision_resolve removes particles during a collision. 
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
    unsigned int rand_seed; ///< seed for random number generator
    /** @} */
    
    
    /**
     * \name Variables related to SimulationArchive 
     * @{
     */
    int    simulationarchive_version;           ///< Version of the SA binary format (1=original/, 2=incremental)
    long   simulationarchive_size_first;        ///< (Deprecated SAV1) Size of the initial binary file in a SA
    long   simulationarchive_size_snapshot;     ///< (Deprecated SAV1) Size of a snapshot in a SA (other than 1st), in bytes
    double simulationarchive_auto_interval;     ///< Current sampling cadence, in code units
    double simulationarchive_auto_walltime;     ///< Current sampling cadence, in wall time
    unsigned long long simulationarchive_auto_step;  ///< Current sampling cadence, in time steps
    double simulationarchive_next;              ///< Next output time (simulation tim or wall time, depending on wether auto_interval or auto_walltime is set)
    unsigned long long simulationarchive_next_step; ///< Next output step (only used if auto_steps is set)
    char*  simulationarchive_filename;          ///< Name of output file
    /** @} */

    /**
     * \name Variables describing the current module selection 
     * @{
     */
    /**
     * @brief Available visualization options
     */
    enum {
        REB_VISUALIZATION_NONE = 0,     ///< No visualization (default if OPENGL compiler flag is turned off)
        REB_VISUALIZATION_OPENGL = 1,   ///< OpenGL visualization (default if OPENGL compiler flag is turned on)
        REB_VISUALIZATION_WEBGL = 2,    ///< WebGL visualization, only usable from Jupyter notebook widget
        } visualization;
    /**
     * @brief Available collision routines
     */
    enum {
        REB_COLLISION_NONE = 0,     ///< Do not search for collisions (default)
        REB_COLLISION_DIRECT = 1,   ///< Direct collision search O(N^2)
        REB_COLLISION_TREE = 2,     ///< Tree based collision search O(N log(N))
        REB_COLLISION_MERCURIUS = 3,///< OBSOLETE, use REB_COLLISION_DIRECT instead
        REB_COLLISION_LINE = 4,     ///< Direct collision search O(N^2), looks for collisions by assuming a linear path over the last timestep
        REB_COLLISION_LINETREE = 5, ///< Tree-based collision search O(N log(N)), looks for collisions by assuming a linear path over the last timestep
        } collision;
    /**
     * @brief Available integrators
     */
    enum {
        REB_INTEGRATOR_IAS15 = 0,    ///< IAS15 integrator, 15th order, non-symplectic (default)
        REB_INTEGRATOR_WHFAST = 1,   ///< WHFast integrator, symplectic, 2nd order, up to 11th order correctors
        REB_INTEGRATOR_SEI = 2,      ///< SEI integrator for shearing sheet simulations, symplectic, needs OMEGA variable
        REB_INTEGRATOR_LEAPFROG = 4, ///< LEAPFROG integrator, simple, 2nd order, symplectic
        // REB_INTEGRATOR_HERMES = 5,   ///< HERMES Integrator, not available anymore. Use MERCURIUS instead.
        REB_INTEGRATOR_NONE = 7,     ///< Do not integrate anything
        REB_INTEGRATOR_JANUS = 8,    ///< Bit-wise reversible JANUS integrator.
        REB_INTEGRATOR_MERCURIUS = 9,///< MERCURIUS integrator 
        REB_INTEGRATOR_SABA = 10,    ///< SABA integrator family (Laskar and Robutel 2001)
        REB_INTEGRATOR_EOS = 11,     ///< Embedded Operator Splitting (EOS) integrator family (Rein 2019)
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
        REB_GRAVITY_MERCURIUS = 4,  ///< Special gravity routine only for MERCURIUS
        REB_GRAVITY_JACOBI = 5,     ///< Special gravity routine which includes the Jacobi terms for WH integrators 
        } gravity;
    /** @} */


    /**
     * \name Integrator structs (the contain integrator specific variables and temporary data structures) 
     * @{
     */
    struct reb_simulation_integrator_sei ri_sei;        ///< The SEI struct 
    struct reb_simulation_integrator_whfast ri_whfast;  ///< The WHFast struct 
    struct reb_simulation_integrator_saba ri_saba;      ///< The SABA struct 
    struct reb_simulation_integrator_ias15 ri_ias15;    ///< The IAS15 struct
    struct reb_simulation_integrator_mercurius ri_mercurius;      ///< The MERCURIUS struct
    struct reb_simulation_integrator_janus ri_janus;    ///< The JANUS struct 
    struct reb_simulation_integrator_eos ri_eos;        ///< The EOS struct 
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
     * @brief This function allows the user to make changes before each timestep.
     */
    void (*pre_timestep_modifications) (struct reb_simulation* const r);
    /**
     * @brief This function allows the user to make changes after each timestep.
     */
    void (*post_timestep_modifications) (struct reb_simulation* const r);
    /**
     * @brief This function is called at the beginning of the simulation and at the end of
     * each timestep.
     */
    void (*heartbeat) (struct reb_simulation* r);
    /**
     * @brief This function is called at the beginning of the simulation and at the end of
     * each timestep.
     */
    void (*display_heartbeat) (struct reb_simulation* r);
    /**
     * @brief Return the coefficient of restitution. By default, it is NULL, assuming a coefficient of 1.
     * @details The velocity of the collision is given to allow for velocity dependent coefficients
     * of restitution.
     */
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v); 
    /**
     * @brief Resolve collision within this function. By default, it is NULL, assuming hard sphere model.
     * @details A return value of 0 indicates that both particles remain in the simulation. A return value of 1 (2) indicates that particle 1 (2) should be removed from the simulation. A return value of 3 indicates that both particles should be removed from the simulation. 
     */
    int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);

    /**
     * @brief Free particle's ap pointer.  Called in reb_remove function.
     */
    void (*free_particle_ap) (struct reb_particle* p);
    /**
     * @brief Called in reb_free_pointers function for any necessary cleanup in external libraries that depend on the simulation structure.
     */
    void (*extras_cleanup) (struct reb_simulation* r);
    /** @} */
    
    /**
     * \name Hooks for external libraries
     * @{
     */
    /**
     * @brief Pointer to connect additional (optional) libraries, e.g., reboundx
     */
    void* extras;
    /** @} */
};

/** @} */

/**
 * @defgroup MainRebFunctions Main functions
 * @details These are the functions that typically need to be called by the user.
 * @{
 */
/**
 * @brief Creates and initialises a REBOUND simulation
 * @details Allocate memory for one reb_simulation structure, initialise all variables 
 * and return the pointer to the reb_simulation structure. This function must be called 
 * before any particles are added.
 */
struct reb_simulation* reb_create_simulation(void);

/**
 * @brief Creates a deep copy of a REBOUND simulation
 * @details All simulation data, including all particle data will be copied. Function pointers
 * need to be set manually afterwards. 
 * Working on the copy will not affect the original simulation.
 */
struct reb_simulation* reb_copy_simulation(struct reb_simulation* r);

/**
 * @brief Compares to REBOUND simulations bit by bit.
 * @return If r1 and r2 are exactly equal to each other then 0 is returned, otherwise 1. The walltime parameter in 
 * the REBOUND struct is ignored in this comparison.
 * @param r1 The first REBOUND simulation to be compared.
 * @param r2 The second REBOUND simulation to be compared.
 * @param output_option Is set to 1, then the output is printed on the screen. If set to 2, only the return value indicates if the simulations are different. 
 */
int reb_diff_simulations(struct reb_simulation* r1, struct reb_simulation* r2, int output_option);

/**
 * @brief Initialize reb_simulation structure.
 *
 * @details Same as reb_create_simulation() but does not allocate memory for structure itself.
 * @param r Structure to be initialized (needs to be allocated externally).
 */
void reb_init_simulation(struct reb_simulation* r);

/**
 * @brief Perform one integration step
 * @details You rarely want to call this function yourself.
 * Use reb_integrate instead.
 * @param r The rebound simulation to be integrated by one step.
 */
void reb_step(struct reb_simulation* const r);

/**
 * @brief Perform a fixed number of integration steps
 * @details You rarely want to call this function yourself.
 * Use reb_integrate instead.
 * @param r The rebound simulation to be integrated.
 * @param N_steps The number of steps to be taken. 
 */
void reb_steps(struct reb_simulation* const r, unsigned int N_steps);

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
 * @param root_nx The number of root boxes in the x direction.
 * @param root_ny The number of root boxes in the y direction.
 * @param root_nz The number of root boxes in the z direction.
 */
void reb_configure_box(struct reb_simulation* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

/**
 * @brief Frees up all space used by a REBOUND simulation and the reb_simulation structure itself.
 * @details The REBOUND simulation is not usable anymore after being passed to this function.
 * @param r The rebound simulation to be freed
 */
void reb_free_simulation(struct reb_simulation* const r);

/**
 * @cond PRIVATE
 */

/**
 * @brief Frees up all space used by a REBOUND simulation, but not the reb_simulation structure itself.
 * @details The REBOUND simulation is not usable anymore after being passed to this function.
 * @param r The rebound simulation to be freed
 */
void reb_free_pointers(struct reb_simulation* const r);
/** @endcond */

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

#ifdef OPENMP
/**
 * @brief Wrapper method to set number of OpenMP threads from python.
 */
void reb_omp_set_num_threads(int num_threads);
#endif // OPENMP

/**
 * @cond PRIVATE
 */

/**
 * @brief Function used to allow binary input.
 */
void reb_reset_temporary_pointers(struct reb_simulation* const r);
/**
 * @brief Function used to allow binary input.
 * @return Returns 1 if one ore more function pointers were not NULL before.
 */
int reb_reset_function_pointers(struct reb_simulation* const r);
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
 * @brief Remove a particle by its hash.
 * @details see examples/removing_particles_from_simulation.
 * @param r The rebound simulation to be considered
 * @param id The hash of the particle to be removed.
 * @param keepSorted If set to 1 keep the particles with indices in the particles array
 * higher than the one with the passed id are all shifted down one position,
 * ensuring the ordering remains. 
 * @return Returns 1 if particle successfully removed,
 * 0 if hash was not found in the particles array.
 */
int reb_remove_by_hash(struct reb_simulation* const r, uint32_t hash, int keepSorted);

/**
 * @brief Get a pointer to a particle by its hash.
 * @details see examples/uniquely_identifying_particles_with_hashes.
 * @param r The rebound simulation to be considered.
 * @param hash The hash of the particle to search for.
 * @return A pointer to the particle if found, NULL otherwise.
*/
struct reb_particle* reb_get_particle_by_hash(struct reb_simulation* const r, uint32_t hash);

/**
 * @brief Run the heartbeat function and check for escaping/colliding particles.
 * @details You rarely want to call this function yourself. It is used internally to 
 * call the function you set to the heartbeat variable in reb_simulation.
 * @param r The rebound simulation to be considered
 */
void reb_run_heartbeat(struct reb_simulation* const r);


/**
 * @brief A force switching function for the MERCURIUS integrator. This function implements 
 * the same polynomial switching function used in MERCURY. 
 */
double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit);           

/**
 * @brief A force switching function for the MERCURIUS integrator. This function implements 
 * an infinitely differentiable switching function. 
 */
double reb_integrator_mercurius_L_infinite(const struct reb_simulation* const r, double d, double dcrit);           

/**
 * @brief Resolve collision by simply halting the integration and setting r->status=REB_EXIT_COLLISION (Default)
 */
int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c);

/**
 * @brief Hardsphere collision resolving routine (default).
 */
int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);

/**
 * @brief Merging collision resolving routine.
 * @details Merges particle with higher index into particle of lower index.
 *          Conserves mass, momentum and volume.  
 */
int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);

/** @} */

/**
 * @defgroup ToolsRebFunctions Helper functions
 * @{
 */
/**
 * @brief Return uniformly distributed random variable in a given range.
 * @param r Rebound simulation. This is required because it contains the state for the random number generator.
 * @param min Minimum value.
 * @param max Maximum value.
 * @return A random variable
 */
double reb_random_uniform(struct reb_simulation* r, double min, double max);

/**
 * @brief Returns a random variable drawn form a powerlaw distribution.
 * @param r Rebound simulation. This is required because it contains the state for the random number generator.
 * @param min Minimum value.
 * @param max Maximum value.
 * @param slope Slope of powerlaw distribution.
 * @return A random variable
 */
double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);

/**
 * @brief Return a random number with normal distribution.
 * @param r Rebound simulation. This is required because it contains the state for the random number generator.
 * @details Algorithm by D.E. Knut, 1997, The Art of Computer Programmin, Addison-Wesley. 
 * @param variance Variance of normal distribution.
 * @return A random variable
 */
double reb_random_normal(struct reb_simulation* r, double variance);

/**
 * @brief Return a random variable drawn form a Rayleigh distribution.  
 * @param r Rebound simulation. This is required because it contains the state for the random number generator.
 * @details Calculated as described on Rayleigh distribution wikipedia page
 * @param sigma Scale parameter.
 * @return A random variable
 */
double reb_random_rayleigh(struct reb_simulation* r, double sigma);

/**
 * @brief Move to the heliocentric frame.
 * @details This function moves all particles to the heliocentric 
 * frame. Note that the simulation will not stay in the heliocentric frame
 * as it is not an inertial frame. Variational particles are not affected
 * by the function.
 * @param r The rebound simulation to be considered
 */
void reb_move_to_hel(struct reb_simulation* const r);

/**
 * @brief Move to center of momentum and center of mass frame.
 * @details This function moves all particles to the center of mass 
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

/**
 * @brief Sets arrays to particle data. 
 * @details This function can be used to quickly access particle data in a serialized form.
 * NULL pointers will not be set.
 * @param r The rebound simulation to be considered
 * @param hash 1D array to hold particle hashes
 * @param mass 1D array to hold particle masses
 * @param radius 1D array to hold particle radii
 * @param xyz 3D array to hold particle positions
 * @param vxvyvz 3D array to hold particle velocities
 * @param xyzvxvyvz 3D array to hold particle positions and velocities
*/
void reb_serialize_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]);

/**
 * @brief Sets particle data to data provided in arrays. 
 * @details This function can be used to quickly set particle data in a serialized form.
 * NULL pointers will not be accessed.
 * @param r The rebound simulation to be considered
 * @param hash 1D array to hold particle hashes
 * @param mass 1D array to hold particle masses
 * @param radius 1D array to hold particle radii
 * @param xyz 3D array to hold particle positions
 * @param vxvyvz 3D array to hold particle velocities
 * @param xyzvxvyvz 3D array to hold particle positions and velocities
*/
void reb_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]);

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
/** @} */

/**
 * @defgroup OutputRebFunctions Output functions
 * List of the built-in output functions for REBOUND
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
 * @details The orbital parameters are calculated in Jacobi coordinates.
 * Particles are assumed to be sorted from the inside out, the central object having index 0. 
 * Each time the function is called N-1 rows are appended to the file with name filename.
 * Each row in the file corresponds to one particle and contains the following columns (tab separated):
 * time, semi-major axis, eccentricity, inclination, Omega (longitude ascending node), 
 * omega (argument of pericenter), lambda (mean longitude), period, f (true anomaly). 
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
void reb_output_binary(struct reb_simulation* r, const char* filename);

/**
 * @brief This function compares two REBOUND simulations and records the difference in a buffer.
 * @details This is used for taking a SimulationArchive Snapshot.
 * @param buf1 The buffer corresponding to the first rebound simulation to be compared
 * @param buf2 The buffer corresponding to the second rebound simulation to be compared
 * @param bufp The buffer which will contain the differences. 
 */
void reb_binary_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep);

/**
 * @brief Same as reb_binary_diff but with more options.
 * @param output_option If set to 0, the differences are written to bufp. If set to 1, printed on the screen. If set to 2, then only the return value indicates any differences.
 * @return 0 is returned if the simulations do not differ (are equal). 1 is return if they differ.
 */
int reb_binary_diff_with_options(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option);

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
void reb_output_binary_positions(struct reb_simulation* r, const char* filename);

/**
 * @brief Append the velocity dispersion of the particles to an ASCII file.
 * @param r The rebound simulation to be considered
 * @param filename Output filename.
 */
void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename);
/** @} */

/**
 * @defgroup SetupRebFunctions Setup functions
 * List of the built-in setup helper functions
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
 * @param err Pointer to error code that will be set by this function. Used for checking why particle was set to nans.
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
 * @details Pal describes a coordinate system for Keplerian Orbits that is analytical (i.e. infinitely differentiable) between spatial coordinates and orbital elements. See http://adsabs.harvard.edu/abs/2009MNRAS.396.1737P
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
 * This can be used to restart a simulation.
 * @param filename Filename to be read.
 * @return Returns a pointer to a REBOUND simulation.
 */
struct reb_simulation* reb_create_simulation_from_binary(char* filename);

/**
 * @brief Enum describing possible errors that might occur during binary file reading.
 */
enum reb_input_binary_messages {
    REB_INPUT_BINARY_WARNING_NONE = 0,
    REB_INPUT_BINARY_ERROR_NOFILE = 1,
    REB_INPUT_BINARY_WARNING_VERSION = 2,
    REB_INPUT_BINARY_WARNING_POINTERS = 4,
    REB_INPUT_BINARY_WARNING_PARTICLES = 8,
    REB_INPUT_BINARY_ERROR_FILENOTOPEN = 16,
    REB_INPUT_BINARY_ERROR_OUTOFRANGE = 32,
    REB_INPUT_BINARY_ERROR_SEEK = 64,
    REB_INPUT_BINARY_WARNING_FIELD_UNKOWN = 128,
    REB_INPUT_BINARY_ERROR_INTEGRATOR = 256,
};

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

/**
 * \name Derivative functions
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
 * @param po The original particle for which the derivative is to be calculated.
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
/** @} */

/**
 * @defgroup ParticleManipFunctions Particle manipulation functions
 * List of reb_particle manipulation functions for REBOUND
 * @{
 */
/**
 * @brief Subtract particle p2 from particle p1 (p1 - p2).
 * @details Subtracts positions, velocities, and mass element by element. 
 * @param p1 First reb_particle (will be modified)
 * @param p2 Second reb_particle to subtract from p1.
 * @returns A new particle with no pointers (not in any simulation etc.) set.
 */
void reb_particle_isub(struct reb_particle* p1, struct reb_particle* p2);

/**
 * @brief Add particle p2 to particle p1.
 * @details Adds positions, velocities, and mass element by element. 
 * @param p1 First reb_particle (will be modified)
 * @param p2 Second reb_particle.
 */
void reb_particle_iadd(struct reb_particle* p1, struct reb_particle* p2);

/**
 * @brief Multiply a particle's members by a constant.
 * @details Multiplies particle's positions, velocities, and mass by a constant.
 * @param p1 reb_particle to modify.
 * @param value Value by which to multiply particle's fields.
 */
void reb_particle_imul(struct reb_particle* p1, double value);

/**
 * @brief Calculate the distance between two particles.
 * @param p1 reb_particle First particle.
 * @param p2 reb_particle Second particle.
 * @param value Distance between p1 and p2.
 */
double reb_particle_distance(struct reb_particle* p1, struct reb_particle* p2);

/** @} */

/**
 * @defgroup SimulationArchiveFunctions Simulation Archive functions
 * Functions for interacting with simulation archives
 * @{
 */

/**
 * @brief Allocates a simulation and sets it to a specific snapshot in a SimulationArchive file.
 */
struct reb_simulation* reb_create_simulation_from_simulationarchive(struct reb_simulationarchive* sa, long snapshot);

/**
 * @brief Equivalent to reb_create_simulation_from_simulationarchive() but also processes warning messages.
 */
void reb_create_simulation_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, long snapshot, enum reb_input_binary_messages* warnings);

/**
 * @brief Opens a SimulationArchive
 * @details This function opens a SimulationArchive file and creates an index to 
 * find snapshots within the file. This may take a few seconds if there are many 
 * snapshots in the file. Note that opening a file on a slow filesystem (for example
 * via a network) might be particularly slow. The file is kept open until 
 * the user calls reb_close_simulationarchive.
 * @param filename The path and filename of the SimulationArchive input file.
 * @return Returns a reb_simulationarchive struct.
 */
struct reb_simulationarchive* reb_open_simulationarchive(const char* filename);

/**
 * @brief Closes a SimulationArchive
 * @details This function closes a SimulationArchive that was previously opened
 * by the reb_open_simulationarchive function. It also destroys the index created
 * and frees the allocated memory. 
 * @param sa The SimulationArchive to be closed.
 */
void reb_close_simulationarchive(struct reb_simulationarchive* sa);

/**
 * @brief Appends a SimulationArchive snapshot to a file
 * @details This function can either be called manually or via one of the convenience methods
 * reb_simulationarchive_automate_interval and reb_simulationarchive_automate_walltime.
 * If the file does not exist, the function outputs a binary file. If a file exists,
 * it appends a SimulationArchive snapshot to that file. 
 * @param r The rebound simulation to be considered.
 * @param filename The path and filename of the SimulationArchive output file.
 */
void reb_simulationarchive_snapshot(struct reb_simulation* r, const char* filename);

/**
 * @brief Automatically create a SimulationArchive Snapshot at regular intervals
 * @param r The rebound simulation to be considered.
 * @param filename The path and filename of the SimulationArchive output file.
 * @param interval The interval between snapshots (in simulation time units).
 */
void reb_simulationarchive_automate_interval(struct reb_simulation* const r, const char* filename, double interval);

/**
 * @brief Automatically create a SimulationArchive Snapshot at regular walltime intervals
 * @param r The rebound simulation to be considered.
 * @param filename The path and filename of the SimulationArchive output file.
 * @param interval The walltime interval between snapshots (in seconds).
 */
void reb_simulationarchive_automate_walltime(struct reb_simulation* const r, const char* filename, double walltime);

/**
 * @brief Automatically create a SimulationArchive Snapshot after a fixed number of timesteps
 * @param r The rebound simulation to be considered.
 * @param filename The path and filename of the SimulationArchive output file.
 * @param unsigned long long The number of timesteps between snapshots.
 */
void reb_simulationarchive_automate_step(struct reb_simulation* const r, const char* filename, unsigned long long step);


/**
 * @cond PRIVATE
 */

/**
 * @brief Frees all the pointers in a SimulationArchive structure
 * @param sa The SimulationArchive to be closed.
 */
void reb_free_simulationarchive_pointers(struct reb_simulationarchive* sa);

/** @endcond */

/** @} */

/**
 * @defgroup TransformationFunctions Coordinate transformations
 * Functions for transforming between various coordinate systems.
 * @{
 */

/**
 * \name Functions to obtain various masses needed in Jacobi coordinates.
 * @{
 */
/**
 * \name From inertial to Jacobi coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_j Particles array in which the Jacobi quantities will be stored.
 * @param p_mass Should be the same particles array as ps for real particles. If passing variational
 * particles in ps, p_mass should be the corresponding array of real particles.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */

void reb_transformations_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j,const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
/** @} */
/**
 * \name From Jacobi to inertial coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_j Particles array in which the Jacobi quantities will be stored.
 * @param p_mass Should be the same particles array as ps for real particles. If passing variational
 * particles in ps, p_mass should be the corresponding array of real particles.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */
void reb_transformations_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
/** @} */
/**
 * \name From inertial to democratic heliocentric coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_h Particles array in which the democratic heliocentric quantities will be stored.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */
void reb_transformations_inertial_to_democraticheliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active);
/** @} */
/**
 * \name From democratic heliocentric to inertial coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_h Particles array in which the democratic heliocentric quantities will be stored.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */
void reb_transformations_democraticheliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_democraticheliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
/** @} */
/**
 * \name From inertial to WHDS coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_h Particles array in which the WHDS quantities will be stored.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */
void reb_transformations_inertial_to_whds_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active);
/** @} */
/**
 * \name From WHDS to inertial coordinates
 * @{
 * @details Different functions allow you to calculate subsets of the positions (pos), velocities(vel),
 *          and accelerations (acc) depending on what's needed.
 * @param ps Particles array with the inertial quantities to convert
 * @param p_h Particles array in which the WHDS quantities will be stored.
 * @param N number of particles in the array.
 * @param N_active number of active particles in the array.
 */
void reb_transformations_whds_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_whds_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
/** @} */
/** @} */

/**
 * @defgroup MiscRebFunctions Miscellaneous functions
 * List of the miscellaneous helper functions for REBOUND
 * @{
 */
/**
 * @brief Calculate the total energy (potential and kinetic).
 * @details Does not work for SEI (shearing sheet simulations). 
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
 * @brief Init the MEGNO particles, enable MEGNO calculation, and specify a seed for the random number generation.
 * @param r The rebound simulation to be considered
 * @param seed The seed to use for the random number generator
 */
void reb_tools_megno_init_seed(struct reb_simulation* const r, unsigned int seed);

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
 * @brief Returns hash for passed string.
 * @param str String key. 
 * @return hash for the passed string.
 */
uint32_t reb_hash(const char* str);

/**
 * @brief Returns a reb_particle structure with fields/hash/ptrs initialized to nan/0/NULL. 
 * @return reb_particle with fields initialized to nan.
 */
struct reb_particle reb_particle_nan(void);

/**
 * @brief Print out an error message, then exit in a semi-nice way.
 */
void reb_exit(const char* const msg);

/**
 * @brief Print or store a warning message, then continue.
 */
void reb_warning(struct reb_simulation* const r, const char* const msg);

/**
 * @brief Print or store an error message, then continue.
 */
void reb_error(struct reb_simulation* const r, const char* const msg);

/**
 * @brief Get the next warning message stored. Used only if save_messages==1.
 * @param r The rebound simulation to be considered
 * @param buf The buffer in which the error message it copied (needs to be at least reb_max_messages_length long).
 * @return Return value is 0 if no messages are present, 1 otherwise.
 */
int reb_get_next_message(struct reb_simulation* const r, char* const buf);
/** @} */
/** @} */

/**
 * @cond PRIVATE
 * Internal functions for calling various integrator steps. Nothing to be changed by the user.
 */

void reb_integrator_whfast_from_inertial(struct reb_simulation* const r);   ///< Internal function to the appropriate WHFast coordinates from inertial
void reb_integrator_whfast_to_inertial(struct reb_simulation* const r); ///< Internal function to move back from particular WHFast coordinates to inertial
void reb_integrator_whfast_reset(struct reb_simulation* r);		///< Internal function used to call a specific integrator
void reb_whfast_interaction_step(struct reb_simulation* const r, const double _dt);///< Internal function
void reb_whfast_jump_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
void reb_whfast_kepler_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
void reb_whfast_com_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
void reb_integrator_ias15_part2(struct reb_simulation* r);              ///< Internal function used to call a specific integrator
void reb_integrator_ias15_reset(struct reb_simulation* r);              ///< Internal function used to call a specific integrator
int reb_integrator_whfast_init(struct reb_simulation* const r);    ///< Internal function to check errors and allocate memory if needed

/** 
 * @brief This function updates the acceleration on all particles. 
 * @details It uses the current position and velocity data in the 
 * (struct reb_particle*) particles structure.
 */
void reb_update_acceleration(struct reb_simulation* r);

/** @endcond */

/**
* @cond PRIVATE
* Related to OpenGL/WebGL visualization. Nothing to be changed by the user.
*/
struct reb_quaternion {
    double x, y, z, w;
};
struct reb_particle_opengl {
    float x,y,z;
    float vx,vy,vz;
    float r;
};
struct reb_orbit_opengl {
    float x,y,z;
    float a, e, f;
    float omega, Omega, inc;
};

struct reb_display_data {
    struct reb_simulation* r;
    struct reb_simulation* r_copy;
    struct reb_particle_opengl* particle_data;
    struct reb_orbit_opengl* orbit_data;
    struct reb_particle* particles_copy;
    struct reb_particle* p_jh_copy;
    unsigned long allocated_N;
    unsigned long allocated_N_whfast;
    unsigned int opengl_enabled;
    double scale;
    double mouse_x;
    double mouse_y;
    double retina;
    pthread_mutex_t mutex;          /**< Mutex to guarantee non-flickering */
    int spheres;                    /**< Switches between point sprite and real spheres. */
    int pause;                      /**< Pauses visualization, but keep simulation running */
    int wire;                       /**< Shows/hides orbit wires. */
    int onscreentext;               /**< Shows/hides onscreen text. */
    int onscreenhelp;               /**< Shows/hides onscreen help. */
    int multisample;                /**< Turn off/on multisampling. */
    int clear;                      /**< Toggles clearing the display on each draw. */
    int ghostboxes;                 /**< Shows/hides ghost boxes. */
    int reference;                  /**< reb_particle used as a reference for centering. */
    unsigned int mouse_action;      
    unsigned int key_mods;      
    struct reb_quaternion view;
    unsigned int simplefont_tex;
    unsigned int simplefont_shader_program;
    unsigned int simplefont_shader_vao;
    unsigned int simplefont_shader_pos_location;
    unsigned int simplefont_shader_ypos_location;
    unsigned int simplefont_shader_scale_location;
    unsigned int simplefont_shader_aspect_location;
    unsigned int simplefont_shader_charval_buffer;
    unsigned int box_shader_program;
    unsigned int box_shader_box_vao;
    unsigned int box_shader_cross_vao;
    unsigned int box_shader_mvp_location;
    unsigned int box_shader_color_location;
    unsigned int point_shader_mvp_location;
    unsigned int point_shader_color_location;
    unsigned int point_shader_program;
    unsigned int point_shader_particle_vao;
    unsigned int sphere_shader_mvp_location;
    unsigned int sphere_shader_program;
    unsigned int sphere_shader_particle_vao;
    unsigned int sphere_shader_vertex_count;
    unsigned int orbit_shader_mvp_location;
    unsigned int orbit_shader_program;
    unsigned int orbit_shader_particle_vao;
    unsigned int orbit_shader_vertex_count;
};
/**
 * @endcond
 */

#ifdef __cplusplus
}
#endif

#undef REBOUND_RESTRICT

#endif
