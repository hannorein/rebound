#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

// WHFast Integrator (Rein & Tamayo 2015)
struct reb_integrator_whfast_hj_state {

#define REB_INTEGRATOR_WHFAST_HJ_KERNEL(X,Y) \
    X(Y, 0, DEFAULT) \
    X(Y, 1, LAZY)
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_WHFAST_HJ_KERNEL)
    } kernel;                                                   // Kernel type. See Rein, Tamayo & Brown 2019 for details.                            
#define REB_INTEGRATOR_WHFAST_HJ_COORDINATES(X,Y) \
    X(Y, 0, JACOBI)                                         /* Jacobi coordinates (default)                   */ \
    X(Y, 1, DEMOCRATICHELIOCENTRIC)                         /* Democratic Heliocentric coordinates            */ \
    X(Y, 2, WHDS)                                           /* WHDS coordinates (Hernandez and Dehnen, 2017)  */ \
    X(Y, 3, BARYCENTRIC)                                    /* Barycentric coordinates                        */ 
    enum {
        REB_GENERATE_ENUM(REB_INTEGRATOR_WHFAST_HJ_COORDINATES)
    } coordinates;                                              // Coordinate system used in Hamiltonian splitting
    unsigned int safe_mode;                                     // 0: Drift Kick Drift scheme (default), 1: combine first and last sub-step.
    unsigned int keep_unsynchronized;                           // 1: continue from unsynchronized state after synchronization 

    // Internal use
    size_t N_allocated;
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
    size_t N_allocated_var;
    struct reb_particle* REB_RESTRICT p_jh_var; // Jacobi coordinates for variational equations
    size_t N_allocated_temp;
    struct reb_particle* REB_RESTRICT p_temp;   // Used for lazy implementer's kernel 
    unsigned int recalculate_coordinates_but_not_synchronized_warning;
};


REB_API void reb_integrator_whfast_hj_from_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates);
REB_API void reb_integrator_whfast_hj_to_inertial(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates);
REB_API int reb_integrator_whfast_hj_init(struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* whfast);    // Used by REBOUNDx
REB_API void reb_integrator_whfast_hj_interaction_step(struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt);
REB_API void reb_integrator_whfast_hj_jump_step(const struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* whfast, const double _dt); // Used by REBOUNDx
REB_API void reb_integrator_whfast_hj_kepler_step(const struct reb_simulation* const r, struct reb_particle* p_jh, int coordinates, const double _dt);
REB_API void reb_integrator_whfast_hj_com_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt);
// Advances one particle forward in a Keplerian orbit for time dt. mu is the gravitational parameter, G*(m+M). r can be NULL unless variational particles are used or warnings are needed.
REB_API void reb_integrator_whfast_hj_kepler_solver(struct reb_particle* const restrict p, double mu, double dt, const struct reb_simulation* const r);

#endif
