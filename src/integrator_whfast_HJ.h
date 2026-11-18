#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

struct reb_integrator_whfast_hj_node {
    size_t primary;     // Index of real particle or barycenter particle
    size_t secondary;   // Index of real particle or barycenter particle
    struct reb_particle barycenter; 
    struct reb_particle q;  // Relative coordinates
};

struct reb_integrator_whfast_hj_state {
    char* tree;
    struct reb_integrator_whfast_hj_node* nodes;
    size_t tree_N;
};

#endif
