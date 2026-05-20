#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

// WHFast Integrator (Rein & Tamayo 2015)
struct reb_integrator_whfast_hj_state {
    // Internal use
    size_t N_allocated;
    struct reb_particle* REB_RESTRICT p_jh;     // Jacobi/heliocentric/WHDS coordinates
};

#endif
