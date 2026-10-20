#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

struct reb_integrator_whfast_hj_node {
    size_t primary;     // Index of real particle or barycenter particle
    size_t secondary;   // Index of real particle or barycenter particle
    size_t barycenter;  // Index of barycenter particle
    struct reb_particle barycenter_p; 
    struct reb_particle q;  // Relative coordinates
};

struct reb_integrator_whfast_hj_state {
    // One entry per binary orbit, in postorder. Child indices below tree_N
    // refer to particles; larger indices refer to nodes[index - tree_N].
    char* tree;
    struct reb_integrator_whfast_hj_node* nodes;
    size_t tree_N;
    size_t N_barycenters;
};

#endif
