#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;
REB_API int reb_integrator_whfast_hj_tree_to_string(struct reb_simulation* const r, char* const buffer, const size_t buffer_size);
// Tree setters assume whfast_hj is selected, a valid hierarchy with each
// particle exactly once, and positive total mass for every binary.
REB_API int reb_integrator_whfast_hj_set_tree(struct reb_simulation* const r, const char* const tree);
REB_API int reb_integrator_whfast_hj_set_binary_plus_particles_tree(struct reb_simulation* const r);

struct reb_integrator_whfast_hj_node {
    size_t primary;
    size_t secondary;
    size_t barycenter;  // Particle index, or tree_N + shared barycenter index.
    double primary_offset;
    double secondary_offset;
    struct reb_particle q;  // Relative coordinates; m is the orbit's total mass.
};

struct reb_integrator_whfast_hj_state {
    // One entry per binary orbit, in postorder. Child indices below tree_N
    // refer to particles; larger indices refer to nodes[index - tree_N].
    struct reb_integrator_whfast_hj_node* nodes;
    struct reb_particle* barycenters; // Only binaries with two nonzero-mass children.
    double* masses;                  // Original particle masses for validation.
    size_t tree_N;
    size_t N_barycenters;

    // Set after a user-supplied tree has been compiled into nodes.
    int given_tree;
};

#endif
