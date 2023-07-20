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

#define REBOUND_RESTRICT restrict

#include <inttypes.h>
#include <stdint.h>
#include <sys/time.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <math.h>
#ifdef AVX512
#include <immintrin.h>
#endif
#ifdef MPI
#include "mpi.h"
#endif // MPI
#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif // GITHASH
#ifndef __GNUC__
#  define  __attribute__(x)  /*DO NOTHING*/
#endif

extern const char* reb_build_str;   ///< Date and time build string.
extern const char* reb_version_str; ///< Version string.
extern const char* reb_githash_str; ///< Current git hash.
extern const char* reb_logo[26];    ///< Logo of rebound. 
extern volatile sig_atomic_t reb_sigint;  ///< Graceful global interrupt handler 

// Forward declarations
struct reb_simulation;
struct reb_display_data;
struct reb_treecell;
struct reb_variational_configuration;

struct reb_particle {
    double x;
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
    double lastcollision;       // Last time the particle had a physical collision.
    struct reb_treecell* c;     // Pointer to the cell the particle is currently in.
    uint32_t hash;              // Hash, can be used to identify particle.
    void* ap;                   // This pointer allows REBOUNDx to add additional properties to the particle.
    struct reb_simulation* sim; // Pointer to the parent simulation.
};

// Generic 3d vector
struct reb_vec3d {
    double x;
    double y;
    double z;
};

// Rotation (Implemented as a quaternion)
struct reb_rotation {
    double ix;
    double iy;
    double iz;
    double r;
};

// Generic pointer with 7 elements, for internal use only (IAS15).
struct reb_dp7 {
    double* REBOUND_RESTRICT p0;
    double* REBOUND_RESTRICT p1;
    double* REBOUND_RESTRICT p2;
    double* REBOUND_RESTRICT p3;
    double* REBOUND_RESTRICT p4;
    double* REBOUND_RESTRICT p5;
    double* REBOUND_RESTRICT p6;
};

struct reb_ghostbox{
    double shiftx;
    double shifty;
    double shiftz;
    double shiftvx;
    double shiftvy;
    double shiftvz;
};

// Integrator structures 
struct reb_simulation_integrator_ias15 {
    double epsilon;
    double min_dt;
    unsigned int epsilon_global;
    unsigned int dt_mode;
   
    // Internal use
    unsigned long iterations_max_exceeded; // Counter how many times the iteration did not converge. 
    int allocatedN;          
    double* REBOUND_RESTRICT at;
    double* REBOUND_RESTRICT x0;
    double* REBOUND_RESTRICT v0;
    double* REBOUND_RESTRICT a0;
    double* REBOUND_RESTRICT csx;
    double* REBOUND_RESTRICT csv;
    double* REBOUND_RESTRICT csa0;

    struct reb_dp7 g;
    struct reb_dp7 b;
    struct reb_dp7 csb;  // Compensated summation storage for b
    struct reb_dp7 e;
    struct reb_dp7 br;   // Used for resetting the b coefficients if a timestep gets rejected
    struct reb_dp7 er;   // Same for e coefficients

    int* map;               // internal map to particles (this is an identity map except when MERCURIUS is used
    int map_allocated_N;    // allocated size for map
};

struct reb_simulation_integrator_mercurius {
    double (*L) (const struct reb_simulation* const r, double d, double dcrit);  
    double hillfac;        
    unsigned int recalculate_coordinates_this_timestep;
    unsigned int recalculate_dcrit_this_timestep;
    unsigned int safe_mode;
   
    // Internal use
    unsigned int is_synchronized;   
    unsigned int mode;              // 0 if WH is operating, 1 if IAS15 is operating.
    unsigned int encounterN;        // Number of particles currently having an encounter
    unsigned int encounterNactive;  // Number of active particles currently having an encounter
    unsigned int tponly_encounter;  // 0 if any encounters are between two massive bodies. 1 if encounters only involve test particles
    unsigned int allocatedN;
    unsigned int allocatedN_additionalforces;
    unsigned int dcrit_allocatedN;  // Current size of dcrit arrays
    double* dcrit;                  // Precalculated switching radii for particles
    struct reb_particle* REBOUND_RESTRICT particles_backup; //  contains coordinates before Kepler step for encounter prediction
    struct reb_particle* REBOUND_RESTRICT particles_backup_additionalforces; // contains coordinates before Kepler step for encounter prediction
    int* encounter_map;             // Map to represent which particles are integrated with ias15
    struct reb_vec3d com_pos;       // Used to keep track of the centre of mass during the timestep
    struct reb_vec3d com_vel;
};

struct reb_simulation_integrator_sei {
    double OMEGA;
    double OMEGAZ;
    // Internal
    double lastdt;      ///< Cached sin(), tan() for this value of dt.
    double sindt;       ///< Cached sin() 
    double tandt;       ///< Cached tan() 
    double sindtz;      ///< Cached sin(), z axis
    double tandtz;      ///< Cached tan(), z axis
};

struct reb_simulation_integrator_saba {
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
    unsigned int safe_mode;
    unsigned int is_synchronized;
    unsigned int keep_unsynchronized;
};

struct reb_simulation_integrator_whfast {
    unsigned int corrector;
    unsigned int corrector2;
    enum {
        REB_WHFAST_KERNEL_DEFAULT = 0,
        REB_WHFAST_KERNEL_MODIFIEDKICK = 1,
        REB_WHFAST_KERNEL_COMPOSITION = 2,
        REB_WHFAST_KERNEL_LAZY = 3,
    }kernel;
    enum {
        REB_WHFAST_COORDINATES_JACOBI = 0,                      ///< Jacobi coordinates (default)
        REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC = 1,      ///< Democratic Heliocentric coordinates
        REB_WHFAST_COORDINATES_WHDS = 2,                        ///< WHDS coordinates (Hernandez and Dehnen, 2017)
        } coordinates;
    unsigned int recalculate_coordinates_this_timestep;
    unsigned int safe_mode;
    unsigned int keep_unsynchronized;
    // Internal 
    struct reb_particle* REBOUND_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    struct reb_particle* REBOUND_RESTRICT p_temp;   // Used for lazy implementer's kernel 
    unsigned int is_synchronized;
    unsigned int allocated_N;
    unsigned int allocated_Ntemp;
    unsigned int timestep_warning;
    unsigned int recalculate_coordinates_but_not_synchronized_warning;
};

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
    double m[8];
    double x[8];
    double y[8];
    double z[8];
    double vx[8];
    double vy[8];
    double vz[8];
#endif // AVX512
};

struct reb_simulation_integrator_whfast512 {
    unsigned int is_synchronized;
    unsigned int keep_unsynchronized;
    unsigned int allocated_N;
    unsigned int gr_potential;
    unsigned int recalculate_constants;
    struct reb_particle_avx512* p_jh;
    struct reb_particle p_jh0;
};

struct reb_ode{ // defines an ODE 
    unsigned int length; // number of components / dimenion
    unsigned int allocatedN;
    unsigned int needs_nbody;
    double* y;      // Current state 
    double* scale;
    double* C;      // Temporary internal array (extrapolation) 
    double** D;     // Temporary internal array (extrapolation) 
    double* y1;     // Temporary internal array (state during the step) 
    double* y0Dot;  // Temporary internal array (derivatives at beginning of step)
    double* yDot;   // Temporary internal array (derivatives)
    double* yTmp;   // Temporary internal array (midpoint method)
    void (*derivatives)(struct reb_ode* const ode, double* const yDot, const double* const y, const double t); // right hand side 
    void (*getscale)(struct reb_ode* const ode, const double* const y0, const double* const y1); // right hand side (optional) 
    void (*pre_timestep)(struct reb_ode* const ode, const double* const y0); // gets called just before the ODE integration (optional)
    void (*post_timestep)(struct reb_ode* const ode, const double* const y0); // gets called just after the ODE integration (optional)
    struct reb_simulation* r; // weak reference to main simulation 
    void* ref;  // pointer to any additional data needed for derivatives
};


struct reb_simulation_integrator_bs {
    struct reb_ode* nbody_ode; //
    int* sequence;      // stepsize sequence
    int* costPerStep;   // overall cost of applying step reduction up to iteration k + 1, in number of calls.
    double* costPerTimeUnit; // cost per unit step.
    double* optimalStep; // optimal steps for each order. 
    double* coeff;    // extrapolation coefficients.
    double eps_abs; // Allowed absolute scalar error.
    double eps_rel; // Allowed relative scalar error.
    double min_dt;
    double max_dt;
    double dt_proposed;
    int firstOrLastStep;
    int previousRejected;
    int targetIter;
    int user_ode_needs_nbody; // Do not set manually. Use needs_nbody in reb_ode instead.
};

typedef struct _StumpfCoefficients
{
  double * __restrict__ c0;
  double * __restrict__ c1;
  double * __restrict__ c2;
  double * __restrict__ c3;
}StumpfCoefficients;

typedef struct UNIVERSAL_VARS
{
  double * __restrict__ t0;
  double * __restrict__ tLast;
  double * uv_csq;
  double * uv_csp;
  double * uv_csv;  
  double * dt;
  double * __restrict__ Q0;
  double * __restrict__ V0;
  double * __restrict__ P0;
  double * __restrict__ Q1;
  double * __restrict__ V1;
  double * __restrict__ P1;  
  double * __restrict__ X;
  double * __restrict__ Q0_norm;
  double * __restrict__ beta;
  double * __restrict__ eta;
  double * __restrict__ zeta;
  double * __restrict__ period;
  double * __restrict__ Xperiod;
  uint32_t stateVectorSize;
  uint32_t controlVectorSize;

  StumpfCoefficients C;

  double mu;  /// G*mCentral

  // Variables for storing classical orbital elements.
  double * e;
  double * a;
  double * h;
  double * h_norm;
  double * peri;
  double * apo;  
}UNIVERSAL_VARS;

typedef struct _controlVars {
    double* __restrict__ p0;
    double* __restrict__ p1;
    double* __restrict__ p2;
    double* __restrict__ p3;
    double* __restrict__ p4;
    double* __restrict__ p5;
    double* __restrict__ p6;
    uint32_t size;
}controlVars;

typedef struct RADAU
{
  // State vectors
  double * dX;
  double * dQ;
  double * dP;
  // Buffers for rectifying into before performing a synchronisation.
  double * Xout;
  double * Qout;
  double * Pout;
  uint32_t * rectifiedArray;
  // Buffer for predictors
  double * predictors;
  // Derivatives at start of the step.
  double * __restrict__ dState0;
  double * __restrict__ ddState0;
  // Intermediate derivatives.
  double * __restrict__ dState;
  double * __restrict__ ddState;
  // Compensated summation arrays for gravity
  double * __restrict__ cs_dState0;
  double * __restrict__ cs_ddState0;
  double * __restrict__ cs_dState;
  double * __restrict__ cs_ddState;
  // Compensated summation arrays for B's.
  controlVars cs_B;
  controlVars cs_B1st;
  // Compensated summation array for the predictor and corrector.
  double * cs_dX;
  double * cs_dq;
  double * cs_dp;
  // Integrator coefficients.
  controlVars G;
  controlVars B;
  controlVars Blast;
  controlVars Blast_1st;
  controlVars G_1st;
  controlVars B_1st;  
  // Variables for performance metrics
  uint64_t fCalls;
  uint64_t rectifications;
  uint32_t convergenceIterations;
  // Iteration convergence variables.
  double * b6_store;
  double * acc_ptr;
}RADAU;

typedef struct DHEM
{
  // Pointers to osculating orbits at a single point in time.
  double * Xosc;
  double * Qosc;
  double * Posc;
  // Osculating orbits for all stages within a step.
  double * XoscStore;
  double ** XoscArr;
  double * XoscPredStore;
  double ** XoscPredArr;  
  // CS variables for the osculating orbits.
  double * Xosc_cs;
  double * XoscStore_cs;
  double ** XoscArr_cs;
  double * Qosc_cs;
  double * Posc_cs;
  // Pointers to osculating orbit derivatives at a single point in time.
  double * Xosc_dot;
  double * Qosc_dot;
  double * Posc_dot;
  // Osculating orbits derivatives for all stages within a step.
  double * Xosc_dotStore;
  double ** Xosc_dotArr;
  // Variables for summing Xosc+dX into.
  double * X;
  double * Q;
  double * P;
  // Mass variables.
  double * __restrict__ m;
  double * __restrict__ m_inv;
  double mTotal;
  // Rectification variables.
  double * rectifyTimeArray;    /// The time at which we need to rectify each body.
  double * rectificationPeriod; /// Elapsed time to trigger a rectification.
}DHEM;

struct reb_simulation_integrator_tes {
    double dq_max;              // value of dq/q that triggers a rectification (backup to recti_per_orbit)
    double recti_per_orbit;     // main method for triggering a rectification 
    double epsilon;             // tolerance parameter
    double orbital_period;      // The lowest initial orbital period.

    // Synchronisation particle storage.
    uint32_t allocated_N;
    struct reb_particle* particles_dh;

    // Vector dimensions variables.
    uint32_t stateVectorLength;		/// Length of the state vector in doubles.
    uint32_t stateVectorSize;		/// Size in bytes of the state vector.
    uint32_t controlVectorLength;	/// Length of the control vector in doubles.
    uint32_t controlVectorSize; 	/// Size in bytes of n * sizeof(double).

    // State storage
    double * mass;					/// Initial particle masses
    double * X_dh;					/// Memory for current state in dh coords.
	double * Q_dh;					/// Current state in dh coords.
	double * P_dh;					/// Current state in dh coords. 
    double COM[3];                  // Centre of mass
    double COM_dot[3];              // Velocity of COM

    // Pointers to various modules comprising TES.
    UNIVERSAL_VARS * uVars;			/// Pointer to the universal variables module
    DHEM * rhs;						/// Pointer to the DHEM rhs
    RADAU * radau;  				/// Pointer to our integrator

    double mStar_last;              /// Mass of the star last step.
    uint32_t warnings;              /// Number of times warning has been shown
};

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

struct reb_simulation_integrator_eos {
    enum REB_EOS_TYPE phi0;
    enum REB_EOS_TYPE phi1;
    unsigned int n;
    unsigned int safe_mode;
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

struct reb_simulation_integrator_janus {
    double scale_pos;
    double scale_vel;
    unsigned int order;
    unsigned int recalculate_integer_coordinates_this_timestep;
    struct reb_particle_int* REBOUND_RESTRICT p_int;
    unsigned int allocated_N;
};

struct reb_collision{
    int p1;
    int p2;
    struct reb_ghostbox gb;
    int ri;
};

// Possible return values of of rebound_integrate
enum REB_STATUS {
    REB_RUNNING_PAUSED = -3,    // Simulation is paused by visualization.
    REB_RUNNING_LAST_STEP = -2, // Current timestep is the last one. Needed to ensure that t=tmax exactly.
    REB_RUNNING = -1,           // Simulation is current running, no error occurred.
    REB_EXIT_SUCCESS = 0,       // Integration finished successfully.
    REB_EXIT_ERROR = 1,         // A generic error occurred and the integration was not successful.
    REB_EXIT_NOPARTICLES = 2,   // The integration ends early because no particles are left in the simulation.
    REB_EXIT_ENCOUNTER = 3,     // The integration ends early because two particles had a close encounter (see exit_min_distance)
    REB_EXIT_ESCAPE = 4,        // The integration ends early because a particle escaped (see exit_max_distance)  
    REB_EXIT_USER = 5,          // User caused exit, simulation did not finish successfully.
    REB_EXIT_SIGINT = 6,        // SIGINT received. Simulation stopped.
    REB_EXIT_COLLISION = 7,     // The integration ends early because two particles collided. 
};

// IDs for content of a binary field. Used to read and write binary files.
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
    REB_BINARY_FIELD_TYPE_IAS15_DTMODE = 401,
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
    REB_BINARY_FIELD_TYPE_MERCURIUS_RECALCULATE_COORD = 123,
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
    REB_BINARY_FIELD_TYPE_RAND_SEED = 154,
    REB_BINARY_FIELD_TYPE_TESTPARTICLEHIDEWARNINGS = 155,
    REB_BINARY_FIELD_TYPE_BS_EPSABS = 156,
    REB_BINARY_FIELD_TYPE_BS_EPSREL = 157,
    REB_BINARY_FIELD_TYPE_BS_MINDT = 158,
    REB_BINARY_FIELD_TYPE_BS_MAXDT = 159,
    REB_BINARY_FIELD_TYPE_BS_FIRSTORLASTSTEP = 160,
    REB_BINARY_FIELD_TYPE_BS_PREVIOUSREJECTED = 161,
    REB_BINARY_FIELD_TYPE_BS_TARGETITER = 162,
    REB_BINARY_FIELD_TYPE_VARRESCALEWARNING = 163,

    REB_BINARY_FIELD_TYPE_TES_DQ_MAX = 300,
    REB_BINARY_FIELD_TYPE_TES_RECTI_PER_ORBIT = 301,
    REB_BINARY_FIELD_TYPE_TES_EPSILON = 302,
    REB_BINARY_FIELD_TYPE_TES_PERIOD = 303,
    REB_BINARY_FIELD_TYPE_TES_ALLOCATED_N = 304,
    REB_BINARY_FIELD_TYPE_TES_PARTICLES_DH = 305,
    REB_BINARY_FIELD_TYPE_TES_SV_LEN = 306,
    REB_BINARY_FIELD_TYPE_TES_SV_SIZE = 307,
    REB_BINARY_FIELD_TYPE_TES_CV_LEN = 308,
    REB_BINARY_FIELD_TYPE_TES_CV_SIZE = 309,
    REB_BINARY_FIELD_TYPE_TES_MASS = 310,
    REB_BINARY_FIELD_TYPE_TES_X_DH = 311,
    REB_BINARY_FIELD_TYPE_TES_COM = 312,
    REB_BINARY_FIELD_TYPE_TES_COM_DOT = 313,
    REB_BINARY_FIELD_TYPE_TES_MASS_STAR_LAST = 314,

    REB_BINARY_FIELD_TYPE_TES_UVARS_SV_SIZE = 320,
    REB_BINARY_FIELD_TYPE_TES_UVARS_CV_SIZE = 321,
    REB_BINARY_FIELD_TYPE_TES_UVARS_T0 = 322,
    REB_BINARY_FIELD_TYPE_TES_UVARS_TLAST = 323,
    REB_BINARY_FIELD_TYPE_TES_UVARS_CSQ = 324,
    REB_BINARY_FIELD_TYPE_TES_UVARS_CSP = 325,
    REB_BINARY_FIELD_TYPE_TES_UVARS_CSV = 326,
    REB_BINARY_FIELD_TYPE_TES_UVARS_Q0 = 327,
    REB_BINARY_FIELD_TYPE_TES_UVARS_V0 = 328,
    REB_BINARY_FIELD_TYPE_TES_UVARS_P0 = 329,
    REB_BINARY_FIELD_TYPE_TES_UVARS_Q1 = 330,
    REB_BINARY_FIELD_TYPE_TES_UVARS_V1 = 331,
    REB_BINARY_FIELD_TYPE_TES_UVARS_P1 = 332,
    REB_BINARY_FIELD_TYPE_TES_UVARS_X = 333,
    REB_BINARY_FIELD_TYPE_TES_UVARS_Q0_NORM = 334,
    REB_BINARY_FIELD_TYPE_TES_UVARS_BETA = 335,
    REB_BINARY_FIELD_TYPE_TES_UVARS_ETA = 336,
    REB_BINARY_FIELD_TYPE_TES_UVARS_ZETA = 337,
    REB_BINARY_FIELD_TYPE_TES_UVARS_PERIOD = 338,
    REB_BINARY_FIELD_TYPE_TES_UVARS_XPERIOD = 339,
    REB_BINARY_FIELD_TYPE_TES_UVARS_STUMPF_C0     = 340,
    REB_BINARY_FIELD_TYPE_TES_UVARS_STUMPF_C1     = 341,
    REB_BINARY_FIELD_TYPE_TES_UVARS_STUMPF_C2     = 342,
    REB_BINARY_FIELD_TYPE_TES_UVARS_STUMPF_C3     = 343,
    REB_BINARY_FIELD_TYPE_TES_UVARS_MU = 344,

    REB_BINARY_FIELD_TYPE_TES_RADAU_DX = 350,
    REB_BINARY_FIELD_TYPE_TES_RADAU_XOUT = 351,
    REB_BINARY_FIELD_TYPE_TES_RADAU_RECTI_ARRAY = 352,
    REB_BINARY_FIELD_TYPE_TES_RADAU_PREDICTORS = 353,
    REB_BINARY_FIELD_TYPE_TES_RADAU_DSTATE0 = 354,
    REB_BINARY_FIELD_TYPE_TES_RADAU_DDSTATE0 = 355,
    REB_BINARY_FIELD_TYPE_TES_RADAU_DSTATE = 356,
    REB_BINARY_FIELD_TYPE_TES_RADAU_DDSTATE = 357,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_DSTATE0 = 358,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_DDSTATE0 = 359,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_DSTATE = 360,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_DDSTATE = 361,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_DX = 362,
    REB_BINARY_FIELD_TYPE_TES_RADAU_FCALLS = 363,
    REB_BINARY_FIELD_TYPE_TES_RADAU_RECTIS = 364,
    REB_BINARY_FIELD_TYPE_TES_RADAU_ITERS = 365,
    REB_BINARY_FIELD_TYPE_TES_RADAU_B6 = 366,    
    REB_BINARY_FIELD_TYPE_TES_RADAU_B = 367,
    REB_BINARY_FIELD_TYPE_TES_RADAU_BLAST = 368,
    REB_BINARY_FIELD_TYPE_TES_RADAU_B_1ST = 369,
    REB_BINARY_FIELD_TYPE_TES_RADAU_BLAST_1ST = 370,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_B = 371,
    REB_BINARY_FIELD_TYPE_TES_RADAU_CS_B_1ST = 372,
    REB_BINARY_FIELD_TYPE_TES_RADAU_G = 373,
    REB_BINARY_FIELD_TYPE_TES_RADAU_G_1ST = 374,

    REB_BINARY_FIELD_TYPE_TES_DHEM_XOSC_STORE = 380,
    REB_BINARY_FIELD_TYPE_TES_DHEM_XOSC_PRED_STORE = 381,
    REB_BINARY_FIELD_TYPE_TES_DHEM_XOSC_CS_STORE = 382,
    REB_BINARY_FIELD_TYPE_TES_DHEM_XOSC_DOT_STORE = 383,
    REB_BINARY_FIELD_TYPE_TES_DHEM_X = 384,
    REB_BINARY_FIELD_TYPE_TES_DHEM_M_INV = 385,
    REB_BINARY_FIELD_TYPE_TES_DHEM_M_TOTAL = 386,
    REB_BINARY_FIELD_TYPE_TES_DHEM_RECTI_TIME = 387,
    REB_BINARY_FIELD_TYPE_TES_DHEM_RECTI_PERIOD = 388,
    
    REB_BINARY_FIELD_TYPE_WHFAST512_KEEPUNSYNC = 390,
    REB_BINARY_FIELD_TYPE_WHFAST512_ISSYNCHRON = 391,
    REB_BINARY_FIELD_TYPE_WHFAST512_GRPOTENTIAL = 392,
    REB_BINARY_FIELD_TYPE_WHFAST512_ALLOCATEDN = 393,
    REB_BINARY_FIELD_TYPE_WHFAST512_PJH = 394,
    REB_BINARY_FIELD_TYPE_WHFAST512_PJH0 = 395,

    REB_BINARY_FIELD_TYPE_HEADER = 1329743186,  // Corresponds to REBO (first characters of header text)
    REB_BINARY_FIELD_TYPE_SABLOB = 9998,        // SA Blob
    REB_BINARY_FIELD_TYPE_END = 9999,
};

// This structure is used to save and load binary files.
struct reb_binary_field {
    uint32_t type;  // enum of REB_BINARY_FIELD_TYPE
    uint64_t size;  // Size in bytes of field (only counting what follows, not the binary field, itself).
};

// Holds a particle's hash and the particle's index in the particles array. Used for particle_lookup_table.
struct reb_hash_pointer_pair{
    uint32_t hash;
    int index;
};

// Main REBOUND Simulation structure
struct reb_simulation {
    double  t;
    double  G;
    double  softening;
    double  dt;
    double  dt_last_done;
    unsigned long long steps_done;
    int     N;
    int     N_var;
    int     var_config_N;
    struct reb_variational_configuration* var_config;   // These configuration structs contain details on variational particles. 
    int     var_rescale_warning;
    int     N_active;
    int     testparticle_type;
    int     testparticle_hidewarnings;
    struct reb_hash_pointer_pair* particle_lookup_table; // Array of pairs that map particles' hashes to their index in the particles array.
    int     hash_ctr;               // Counter for number of assigned hashes to assign unique values.
    int     N_lookup;               // Number of entries in the particle lookup table.
    int     allocatedN_lookup;      // Number of lookup table entries allocated.
    int     allocatedN;             // Current maximum space allocated in the particles array on this node. 
    struct reb_particle* particles;
    struct reb_vec3d* gravity_cs;   // Containing the information for compensated gravity summation 
    int     gravity_cs_allocatedN;
    struct reb_treecell** tree_root;// Pointer to the roots of the trees. 
    int     tree_needs_update;      // Flag to force a tree update (after boundary check)
    double opening_angle2;
    enum REB_STATUS status;
    int     exact_finish_time;

    unsigned int force_is_velocity_dependent;
    unsigned int gravity_ignore_terms;
    double output_timing_last;      // Time when reb_output_timing() was called the last time. 
    unsigned long display_clock;    // Display clock, internal variable for timing refreshs.
    int save_messages;              // Set to 1 to ignore messages (used in python interface).
    char** messages;                // Array of strings containing last messages (only used if save_messages==1). 
    double exit_max_distance;
    double exit_min_distance;
    double usleep;
    struct reb_display_data* display_data; // Datastructure stores visualization related data. Does not have to be modified by the user. 
    int track_energy_offset;
    double energy_offset;
    double walltime;
    uint32_t python_unit_l;         // Only used for when working with units in python.
    uint32_t python_unit_m;         // Only used for when working with units in python.
    uint32_t python_unit_t;         // Only used for when working with units in python.
    
    // Ghost boxes 
    struct  reb_vec3d boxsize;      // Size of the entire box, root_x*boxsize. 
    double  boxsize_max;            // Maximum size of the entire box in any direction. Set in box_init().
    double  root_size;              // Size of a root box. 
    int     root_n;                 // Total number of root boxes in all directions, root_nx*root_ny*root_nz. Default: 1. Set in box_init().
    int     root_nx;                // Number of ghost boxes in x direction. Do not change manually.
    int     root_ny;
    int     root_nz;
    int     nghostx;
    int     nghosty;
    int     nghostz;

#ifdef MPI
    int    mpi_id;                              // Unique id of this node (starting at 0). Used for MPI only.
    int    mpi_num;                             // Number of MPI nodes. Used for MPI only.
    struct reb_particle** particles_send;       // Send buffer for particles. There is one buffer per node. 
    int*   particles_send_N;                    // Current length of particle send buffer. 
    int*   particles_send_Nmax;                 // Maximal length of particle send beffer before realloc() is needed. 
    struct reb_particle** particles_recv;       // Receive buffer for particles. There is one buffer per node. 
    int*   particles_recv_N;                    // Current length of particle receive buffer. 
    int*   particles_recv_Nmax;                 // Maximal length of particle receive beffer before realloc() is needed. */

    struct reb_treecell** tree_essential_send;  // Send buffer for cells. There is one buffer per node. 
    int*   tree_essential_send_N;               // Current length of cell send buffer. 
    int*   tree_essential_send_Nmax;            // Maximal length of cell send beffer before realloc() is needed. 
    struct reb_treecell** tree_essential_recv;  // Receive buffer for cells. There is one buffer per node. 
    int*   tree_essential_recv_N;               // Current length of cell receive buffer. 
    int*   tree_essential_recv_Nmax;            // Maximal length of cell receive beffer before realloc() is needed. 
#endif // MPI

    int collision_resolve_keep_sorted;
    struct reb_collision* collisions;       ///< Array of all collisions. 
    int collisions_allocatedN;
    double minimum_collision_velocity;
    double collisions_plog;
    double max_radius[2];               // Two largest particle radii, set automatically, needed for collision search.
    long collisions_Nlog;
    
    // MEGNO
    int calculate_megno;    // Do not change manually. Internal flag that determines if megno is calculated (default=0, but megno_init() sets it to the index of variational particles used for megno)
    double megno_Ys;        // Running megno sum (internal use)
    double megno_Yss;       // Running megno sum (internal use)
    double megno_cov_Yt;    // covariance of MEGNO Y and t
    double megno_var_t;     // variance of t 
    double megno_mean_t;    // mean of t
    double megno_mean_Y;    // mean of MEGNO Y
    long   megno_n;         // number of covariance updates
    unsigned int rand_seed; // seed for random number generator
    
     // SimulationArchive 
    int    simulationarchive_version;               // Version of the SA binary format (1=original/, 2=incremental)
    long   simulationarchive_size_first;            // (Deprecated SAV1) Size of the initial binary file in a SA
    long   simulationarchive_size_snapshot;         // (Deprecated SAV1) Size of a snapshot in a SA (other than 1st), in bytes
    double simulationarchive_auto_interval;         // Current sampling cadence, in code units
    double simulationarchive_auto_walltime;         // Current sampling cadence, in wall time
    unsigned long long simulationarchive_auto_step; // Current sampling cadence, in time steps
    double simulationarchive_next;                  // Next output time (simulation tim or wall time, depending on wether auto_interval or auto_walltime is set)
    unsigned long long simulationarchive_next_step; // Next output step (only used if auto_steps is set)
    char*  simulationarchive_filename;              // Name of output file

    // Modules
    enum {
        REB_VISUALIZATION_NONE = 0,     // No visualization (default if OPENGL compiler flag is turned off)
        REB_VISUALIZATION_OPENGL = 1,   // OpenGL visualization (default if OPENGL compiler flag is turned on)
        REB_VISUALIZATION_WEBGL = 2,    // WebGL visualization, only usable from Jupyter notebook widget
        } visualization;
    enum {
        REB_COLLISION_NONE = 0,     // Do not search for collisions (default)
        REB_COLLISION_DIRECT = 1,   // Direct collision search O(N^2)
        REB_COLLISION_TREE = 2,     // Tree based collision search O(N log(N))
        REB_COLLISION_LINE = 4,     // Direct collision search O(N^2), looks for collisions by assuming a linear path over the last timestep
        REB_COLLISION_LINETREE = 5, // Tree-based collision search O(N log(N)), looks for collisions by assuming a linear path over the last timestep
        } collision;
    enum {
        REB_INTEGRATOR_IAS15 = 0,    // IAS15 integrator, 15th order, non-symplectic (default)
        REB_INTEGRATOR_WHFAST = 1,   // WHFast integrator, symplectic, 2nd order, up to 11th order correctors
        REB_INTEGRATOR_SEI = 2,      // SEI integrator for shearing sheet simulations, symplectic, needs OMEGA variable
        REB_INTEGRATOR_LEAPFROG = 4, // LEAPFROG integrator, simple, 2nd order, symplectic
        REB_INTEGRATOR_NONE = 7,     // Do not integrate anything
        REB_INTEGRATOR_JANUS = 8,    // Bit-wise reversible JANUS integrator.
        REB_INTEGRATOR_MERCURIUS = 9,// MERCURIUS integrator 
        REB_INTEGRATOR_SABA = 10,    // SABA integrator family (Laskar and Robutel 2001)
        REB_INTEGRATOR_EOS = 11,     // Embedded Operator Splitting (EOS) integrator family (Rein 2019)
        REB_INTEGRATOR_BS = 12,      // Gragg-Bulirsch-Stoer 
        REB_INTEGRATOR_TES = 20,     // Terrestrial Exoplanet Simulator (TES) 
        REB_INTEGRATOR_WHFAST512 = 21,   // WHFast integrator, optimized for AVX512
        } integrator;
    enum {
        REB_BOUNDARY_NONE = 0,      // Do not check for anything (default)
        REB_BOUNDARY_OPEN = 1,      // Open boundary conditions. Removes particles if they leave the box 
        REB_BOUNDARY_PERIODIC = 2,  // Periodic boundary conditions
        REB_BOUNDARY_SHEAR = 3,     // Shear periodic boundary conditions, needs OMEGA variable
        } boundary;
    enum {
        REB_GRAVITY_NONE = 0,       // Do not calculate graviational forces
        REB_GRAVITY_BASIC = 1,      // Basic O(N^2) direct summation algorithm, choose this for shearing sheet and periodic boundary conditions
        REB_GRAVITY_COMPENSATED = 2,// Direct summation algorithm O(N^2) but with compensated summation, slightly slower than BASIC but more accurate
        REB_GRAVITY_TREE = 3,       // Use the tree to calculate gravity, O(N log(N)), set opening_angle2 to adjust accuracy.
        REB_GRAVITY_MERCURIUS = 4,  // Special gravity routine only for MERCURIUS
        REB_GRAVITY_JACOBI = 5,     // Special gravity routine which includes the Jacobi terms for WH integrators 
        } gravity;

    // Integrators
    struct reb_simulation_integrator_sei ri_sei;            // The SEI struct 
    struct reb_simulation_integrator_whfast ri_whfast;      // The WHFast struct 
    struct reb_simulation_integrator_whfast512 ri_whfast512;      // The WHFast512 struct 
    struct reb_simulation_integrator_saba ri_saba;          // The SABA struct 
    struct reb_simulation_integrator_ias15 ri_ias15;        // The IAS15 struct
    struct reb_simulation_integrator_mercurius ri_mercurius;// The MERCURIUS struct
    struct reb_simulation_integrator_janus ri_janus;        // The JANUS struct 
    struct reb_simulation_integrator_eos ri_eos;            // The EOS struct 
    struct reb_simulation_integrator_bs ri_bs;              // The BS struct
    struct reb_simulation_integrator_tes ri_tes;            // TES struct

    // ODEs
    struct reb_ode** odes;  // all ode sets (includes nbody if BS set as integrator)
    int odes_N;            // number of ode sets
    int odes_allocatedN;   // number of ode sets allocated
    int ode_warnings;

     // Callback functions
    void (*additional_forces) (struct reb_simulation* const r);
    void (*pre_timestep_modifications) (struct reb_simulation* const r);    // used by REBOUNDx
    void (*post_timestep_modifications) (struct reb_simulation* const r);   // used by REBOUNDx
    void (*heartbeat) (struct reb_simulation* r);
    void (*display_heartbeat) (struct reb_simulation* r);
    double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v); 
    int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision);
    void (*free_particle_ap) (struct reb_particle* p);   // used by REBOUNDx 
    void (*extras_cleanup) (struct reb_simulation* r);
    void* extras; // Pointer to connect additional (optional) libraries, e.g., reboundx
};


// Structure representing a Keplerian orbit.
struct reb_orbit {
    double d;        // Radial distance from central object
    double v;        // velocity relative to central object's velocity
    double h;        // Specific angular momentum
    double P;        // Orbital period
    double n;        // Mean motion
    double a;        // Semi-major axis
    double e;        // Eccentricity
    double inc;      // Inclination
    double Omega;    // Longitude of ascending node
    double omega;    // Argument of pericenter
    double pomega;   // Longitude of pericenter
    double f;        // True anomaly
    double M;        // Mean anomaly
    double l;        // Mean Longitude
    double theta;    // True Longitude
    double T;        // Time of pericenter passage
    double rhill;    // Circular Hill radius 
    double pal_h;    // Cartesian component of the eccentricity, h = e*sin(pomega) 
    double pal_k;    // Cartesian component of the eccentricity, k = e*cos(pomega) 
    double pal_ix;   // Cartesian component of the inclination, ix = 2*sin(i/2)*cos(Omega)
    double pal_iy;    // Cartesian component of the inclination, ix = 2*sin(i/2)*sin(Omega)
    struct reb_vec3d hvec;  // specific angular momentum vector
    struct reb_vec3d evec;  // eccentricity vector (mag=ecc, points toward peri)
};


// Simulation life cycle
struct reb_simulation* reb_create_simulation(void);     // allocates memory, then calls reb_init_simulation
void reb_init_simulation(struct reb_simulation* r);    
void reb_free_simulation(struct reb_simulation* const r);
struct reb_simulation* reb_copy_simulation(struct reb_simulation* r);
void reb_free_pointers(struct reb_simulation* const r);
void reb_reset_temporary_pointers(struct reb_simulation* const r);
int reb_reset_function_pointers(struct reb_simulation* const r); // Returns 1 if one ore more function pointers were not NULL before.
// Configure the boundary/root box
void reb_configure_box(struct reb_simulation* const r, const double boxsize, const int root_nx, const int root_ny, const int root_nz);

// Messages and control functions
void reb_exit(const char* const msg); // Print out an error message, then exit in a semi-nice way.
void reb_warning(struct reb_simulation* const r, const char* const msg);   // Print or store a warning message, then continue.
void reb_error(struct reb_simulation* const r, const char* const msg);     // Print or store an error message, then continue.
int reb_get_next_message(struct reb_simulation* const r, char* const buf); // Get the next stored warning message. Used only if save_messages==1. Return value is 0 if no messages are present, 1 otherwise.

// Timestepping
void reb_step(struct reb_simulation* const r);
void reb_steps(struct reb_simulation* const r, unsigned int N_steps);
enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax);
void reb_integrator_synchronize(struct reb_simulation* r);
void reb_integrator_reset(struct reb_simulation* r);
void reb_update_acceleration(struct reb_simulation* r);
void reb_stop(struct reb_simulation* const r); // Stop current integration

// Compare simulations
// If r1 and r2 are exactly equal to each other then 0 is returned, otherwise 1. Walltime is ignored.
// If output_option=1, then output is printed on the screen. If 2, only return value os given. 
int reb_diff_simulations(struct reb_simulation* r1, struct reb_simulation* r2, int output_option);

// Mercurius switching functions
double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit);
double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit);
double reb_integrator_mercurius_L_C4(const struct reb_simulation* const r, double d, double dcrit);
double reb_integrator_mercurius_L_C5(const struct reb_simulation* const r, double d, double dcrit);

// Collision resolve functions
int reb_collision_resolve_halt(struct reb_simulation* const r, struct reb_collision c);
int reb_collision_resolve_hardsphere(struct reb_simulation* const r, struct reb_collision c);
int reb_collision_resolve_merge(struct reb_simulation* const r, struct reb_collision c);

// Random sampling
double reb_random_uniform(struct reb_simulation* r, double min, double max);
double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
double reb_random_normal(struct reb_simulation* r, double variance);
double reb_random_rayleigh(struct reb_simulation* r, double sigma);

// Serialization functions.
void reb_serialize_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // NULL pointers will not be set.
void reb_set_serialized_particle_data(struct reb_simulation* r, uint32_t* hash, double* m, double* radius, double (*xyz)[3], double (*vxvyvz)[3], double (*xyzvxvyvz)[6]); // Null pointers will be ignored.

// Output functions
int reb_output_check(struct reb_simulation* r, double interval);
void reb_output_timing(struct reb_simulation* r, const double tmax);
void reb_output_orbits(struct reb_simulation* r, char* filename);
void reb_output_binary(struct reb_simulation* r, const char* filename);
void reb_output_ascii(struct reb_simulation* r, char* filename);
void reb_output_binary_positions(struct reb_simulation* r, const char* filename);
void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename);

// Compares two simulations, stores difference in buffer.
void reb_binary_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep); 
// Same as reb_binary_diff, but with options.
// output_option If set to 0, the differences are written to bufp. If set to 1, printed on the screen. If set to 2, then only the return value indicates any differences.
// returns 0 is returned if the simulations do not differ (are equal). 1 is return if they differ.
int reb_binary_diff_with_options(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option);

// Input functions
struct reb_simulation* reb_create_simulation_from_binary(char* filename);

// Possible errors that might occur during binary file reading.
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
    REB_INPUT_BINARY_WARNING_CORRUPTFILE = 512,
};

// ODE functions
struct reb_ode* reb_create_ode(struct reb_simulation* r, unsigned int length);
void reb_free_ode(struct reb_ode* ode);

// Miscellaneous functions
uint32_t reb_hash(const char* str);
double reb_tools_mod2pi(double f);
double reb_tools_M_to_f(double e, double M); // True anomaly for a given eccentricity and mean anomaly
double reb_tools_E_to_f(double e, double M); // True anomaly for a given eccentricity and eccentric anomaly
double reb_tools_M_to_E(double e, double M); // Eccentric anomaly for a given eccentricity and mean anomaly
void reb_tools_init_plummer(struct reb_simulation* r, int _N, double M, double R); // This function sets up a Plummer sphere, N=number of particles, M=total mass, R=characteristic radius
void reb_run_heartbeat(struct reb_simulation* const r);  // used internally

// transformations to/from vec3d
struct reb_vec3d reb_tools_spherical_to_xyz(const double mag, const double theta, const double phi);
void reb_tools_xyz_to_spherical(struct reb_vec3d const xyz, double* mag, double* theta, double* phi);


// Functions to add and initialize particles
struct reb_particle reb_particle_nan(void); // Returns a reb_particle structure with fields/hash/ptrs initialized to nan/0/NULL. 
void reb_add(struct reb_simulation* const r, struct reb_particle pt);
void reb_add_fmt(struct reb_simulation* r, const char* fmt, ...);
struct reb_particle reb_particle_new(struct reb_simulation* r, const char* fmt, ...);    // Same as reb_add_fmt() but returns the particle instead of adding it to the simualtion.
struct reb_particle reb_tools_orbit_to_particle_err(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f, int* err);
struct reb_particle reb_tools_orbit_to_particle(double G, struct reb_particle primary, double m, double a, double e, double i, double Omega, double omega, double f);
struct reb_particle reb_tools_pal_to_particle(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy);

// Functions to access and remove particles
void reb_remove_all(struct reb_simulation* const r);
int reb_remove(struct reb_simulation* const r, int index, int keepSorted);
int reb_remove_by_hash(struct reb_simulation* const r, uint32_t hash, int keepSorted);
struct reb_particle* reb_get_particle_by_hash(struct reb_simulation* const r, uint32_t hash);
struct reb_particle reb_get_remote_particle_by_hash(struct reb_simulation* const r, uint32_t hash);
int reb_get_particle_index(struct reb_particle* p); // Returns a particle's index in the simulation it's in. Needs to be in the simulation its sim pointer is pointing to. Otherwise -1 returned.
struct reb_particle reb_get_jacobi_com(struct reb_particle* p); // Returns the Jacobi center of mass for a given particle. Used by python. Particle needs to be in a simulation.

// Orbit calculation
struct reb_orbit reb_tools_particle_to_orbit_err(double G, struct reb_particle p, struct reb_particle primary, int* err);
struct reb_orbit reb_tools_particle_to_orbit(double G, struct reb_particle p, struct reb_particle primary);

// Chaos indicators
void reb_tools_megno_init(struct reb_simulation* const r);
void reb_tools_megno_init_seed(struct reb_simulation* const r, unsigned int seed);
double reb_tools_calculate_megno(struct reb_simulation* r);
double reb_tools_calculate_lyapunov(struct reb_simulation* r);

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
    double lrescale;             // Accumulates the logarithm of rescalings
};

// Add and initialize a set of first order variational particles
// If testparticle is >= 0, then only one variational particle (the test particle) will be added.
// If testparticle is -1, one variational particle for each real particle will be added.
// Returns the index of the first variational particle added
int reb_add_var_1st_order(struct reb_simulation* const r, int testparticle);

// Add and initialize a set of second order variational particles
// Note: a set of second order variational particles requires two sets of first order variational equations.
// If testparticle is >= 0, then only one variational particle (the test particle) will be added.
// If testparticle is -1, one variational particle for each real particle will be added.
// index_1st_order_a is the index of the corresponding first variational particles.
// index_1st_order_b is the index of the corresponding first variational particles.
// Returns the index of the first variational particle added
int reb_add_var_2nd_order(struct reb_simulation* const r, int testparticle, int index_1st_order_a, int index_1st_order_b);

// Rescale all sets of variational particles if their size gets too large (>1e100).
// This can prevent an overflow in floating point numbers. The logarithm of the rescaling
// factor is stored in the reb_variational_configuration's lrescale variable. 
// This function is called automatically every timestep. To avoid automatic rescaling,
// set the reb_variational_configuration's lrescale variable to -1.
// For this function to work, the positions and velocities needs to be synchronized. 
// A warning is presented if the integrator is not synchronized. 
void reb_var_rescale(struct reb_simulation* const r);

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

// Functions to operate on particles
void reb_particle_isub(struct reb_particle* p1, struct reb_particle* p2);
void reb_particle_iadd(struct reb_particle* p1, struct reb_particle* p2);
void reb_particle_imul(struct reb_particle* p1, double value);
double reb_particle_distance(struct reb_particle* p1, struct reb_particle* p2);

// Functions to operate on simulations
void reb_simulation_imul(struct reb_simulation* r, double scalar_pos, double scalar_vel);
int reb_simulation_iadd(struct reb_simulation* r, struct reb_simulation* r2);
int reb_simulation_isub(struct reb_simulation* r, struct reb_simulation* r2);
void reb_move_to_hel(struct reb_simulation* const r);
void reb_move_to_com(struct reb_simulation* const r);

// Diangnostic functions
double reb_tools_energy(const struct reb_simulation* const r);
struct reb_vec3d reb_tools_angular_momentum(const struct reb_simulation* const r);
struct reb_particle reb_get_com(struct reb_simulation* r);
struct reb_particle reb_get_com_of_pair(struct reb_particle p1, struct reb_particle p2);
struct reb_particle reb_get_com_range(struct reb_simulation* r, int first, int last);


// Simulation Archive
struct reb_simulationarchive_blob {  // Used in the binary file to identify data blobs
    int32_t index;                   // Index of previous blob (binary file is 0, first blob is 1)
    int32_t offset_prev;             // Offset to beginning of previous blob (size of previous blob).
    int32_t offset_next;             // Offset to end of following blob (size of following blob).
};
struct reb_simulationarchive_blob16 {  // For backwards compatability only. Will be removed in a future release. 
    int32_t index;
    int16_t offset_prev;
    int16_t offset_next;
};

struct reb_simulationarchive{
    FILE* inf;                   // File pointer (will be kept open)
    char* filename;              // Filename of open file
    int version;                 // SimulationArchive version
    long size_first;             // Size of first snapshot (only used for version 1)
    long size_snapshot;          // Size of snapshot (only used for version 1)
    double auto_interval;        // Interval setting used to create SA (if used)
    double auto_walltime;        // Walltime setting used to create SA (if used)
    unsigned long long auto_step;// Steps in-between SA snapshots (if used)
    long nblobs;                 // Total number of snapshots (including initial binary)
    uint64_t* offset64;            // Index of offsets in file (length nblobs)
    double* t;                   // Index of simulation times in file (length nblobs)
};
struct reb_simulation* reb_create_simulation_from_simulationarchive(struct reb_simulationarchive* sa, long snapshot);
void reb_create_simulation_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, long snapshot, enum reb_input_binary_messages* warnings);
struct reb_simulationarchive* reb_open_simulationarchive(const char* filename);
void reb_close_simulationarchive(struct reb_simulationarchive* sa);
void reb_simulationarchive_snapshot(struct reb_simulation* r, const char* filename);
void reb_simulationarchive_automate_interval(struct reb_simulation* const r, const char* filename, double interval);
void reb_simulationarchive_automate_walltime(struct reb_simulation* const r, const char* filename, double walltime);
void reb_simulationarchive_automate_step(struct reb_simulation* const r, const char* filename, unsigned long long step);
void reb_free_simulationarchive_pointers(struct reb_simulationarchive* sa);


// Functions to convert between coordinate systems

// Jacobi
// p_mass: Should be the same particles array as ps for real particles. If passing variational
//         particles in ps, p_mass should be the corresponding array of real particles.
void reb_transformations_inertial_to_jacobi_posvel(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_inertial_to_jacobi_posvelacc(const struct reb_particle* const particles, struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_inertial_to_jacobi_acc(const struct reb_particle* const particles, struct reb_particle* const p_j,const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_jacobi_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_jacobi_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);
void reb_transformations_jacobi_to_inertial_acc(struct reb_particle* const particles, const struct reb_particle* const p_j, const struct reb_particle* const p_mass, const unsigned int N, const int N_active);

// Democratic heliocentric coordinates
void reb_transformations_inertial_to_democraticheliocentric_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_democraticheliocentric_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_democraticheliocentric_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);

// WHDS
void reb_transformations_inertial_to_whds_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_whds_to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);
void reb_transformations_whds_to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const unsigned int N, const int N_active);

// Rotations
struct reb_rotation reb_rotation_inverse(const struct reb_rotation q);
struct reb_rotation reb_rotation_mul(const struct reb_rotation p, const struct reb_rotation q);

struct reb_rotation reb_rotation_identity();
struct reb_rotation reb_rotation_normalize(const struct reb_rotation q);
struct reb_rotation reb_rotation_conjugate(const struct reb_rotation q);
struct reb_rotation reb_rotation_init_angle_axis(const double angle, struct reb_vec3d axis);
struct reb_rotation reb_rotation_init_from_to(struct reb_vec3d from, struct reb_vec3d to);
struct reb_rotation reb_rotation_init_orbit(const double Omega, const double inc, const double omega);
struct reb_rotation reb_rotation_init_to_new_axes(struct reb_vec3d newz, struct reb_vec3d newx);

struct reb_vec3d reb_vec3d_mul(const struct reb_vec3d v, const double s);
struct reb_vec3d reb_vec3d_add(const struct reb_vec3d v, const struct reb_vec3d w);
double reb_vec3d_length_squared(const struct reb_vec3d v);
double reb_vec3d_dot(const struct reb_vec3d a, const struct reb_vec3d b);
struct reb_vec3d reb_vec3d_cross(const struct reb_vec3d a, const struct reb_vec3d b);
struct reb_vec3d reb_vec3d_normalize(const struct reb_vec3d v);
struct reb_vec3d reb_vec3d_rotate(struct reb_vec3d v, const struct reb_rotation q);
void reb_vec3d_irotate(struct reb_vec3d *v, const struct reb_rotation q);
void reb_particle_irotate(struct reb_particle* p, const struct reb_rotation q);
void reb_simulation_irotate(struct reb_simulation* const sim, const struct reb_rotation q);

void reb_rotation_to_orbital(struct reb_rotation q, double* Omega, double* inc, double* omega);

#ifdef MPI
void reb_mpi_init(struct reb_simulation* const r);
void reb_mpi_finalize(struct reb_simulation* const r);
#endif // MPI

#ifdef OPENMP
// Wrapper method to set number of OpenMP threads from python.
void reb_omp_set_num_threads(int num_threads);
#endif // OPENMP

// The following stuctures are related to OpenGL/WebGL visualization. Nothing to be changed by the user.
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
    pthread_mutex_t mutex;          // Mutex to guarantee non-flickering
    int spheres;                    // Switches between point sprite and real spheres.
    int pause;                      // Pauses visualization, but keep simulation running
    int wire;                       // Shows/hides orbit wires.
    int onscreentext;               // Shows/hides onscreen text.
    int onscreenhelp;               // Shows/hides onscreen help.
    int multisample;                // Turn off/on multisampling.
    int clear;                      // Toggles clearing the display on each draw.
    int ghostboxes;                 // Shows/hides ghost boxes.
    int reference;                  // reb_particle used as a reference for centering.
    unsigned int mouse_action;      
    unsigned int key_mods;      
    struct reb_rotation view;
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


// Temporary. Function declarations needed by REBOUNDx 
void reb_integrator_ias15_reset(struct reb_simulation* r);         ///< Internal function used to call a specific integrator
void reb_integrator_ias15_part2(struct reb_simulation* r);         ///< Internal function used to call a specific integrator
void reb_integrator_whfast_from_inertial(struct reb_simulation* const r);   ///< Internal function to the appropriate WHFast coordinates from inertial
void reb_integrator_whfast_to_inertial(struct reb_simulation* const r); ///< Internal function to move back from particular WHFast coordinates to inertial
void reb_integrator_whfast_reset(struct reb_simulation* r);		///< Internal function used to call a specific integrator
int reb_integrator_whfast_init(struct reb_simulation* const r);    ///< Internal function to check errors and allocate memory if needed
void reb_whfast_interaction_step(struct reb_simulation* const r, const double _dt);///< Internal function
void reb_whfast_jump_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
void reb_whfast_kepler_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
void reb_whfast_com_step(const struct reb_simulation* const r, const double _dt); ///< Internal function
#endif // _MAIN_H
