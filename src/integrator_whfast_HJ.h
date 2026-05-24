#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

struct hj_node
{
    struct reb_particle barycenter_particle;

    struct hj_node *children;
    int N_children;
    int N_allocated;

    int particle_index; // -1 for wrapper, >= 0 for leaf
};

struct reb_integrator_whfast_hj_state {
    // Internal use
    struct hj_node root;
};

#endif
