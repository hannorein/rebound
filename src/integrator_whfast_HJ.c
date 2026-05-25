
#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include <math.h>
#include "transformations.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast_HJ.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))   ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a))   ///< Returns the minimum of a and b

void* reb_integrator_whfast_hj_create();
void reb_integrator_whfast_hj_free(void* state);
void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state);
static void reb_integrator_whfast_hj_node_free(struct hj_node* const node);

const struct reb_integrator reb_integrator_whfast_hj = {
    .step = reb_integrator_whfast_hj_step,
    .create = reb_integrator_whfast_hj_create,
    .free = reb_integrator_whfast_hj_free,
};

// create and free
void* reb_integrator_whfast_hj_create(){
    // Allocate memory and set default parameters.
    struct reb_integrator_whfast_hj_state* whfast = calloc(sizeof(struct reb_integrator_whfast_hj_state),1);
    return whfast;
}

static void reb_integrator_whfast_hj_node_free(struct hj_node *const node)
{
    for (int i = 0; i < node->N_children; i++)
    {
        reb_integrator_whfast_hj_node_free(&node->children[i]);
    }
    free(node->children);
    *node = (struct hj_node){0};
    node->particle_index = -1;
}

void reb_integrator_whfast_hj_free(void* p){
    struct reb_integrator_whfast_hj_state* whfast = p;
    reb_integrator_whfast_hj_node_free(&whfast->root);
    free(whfast);
}


// helper function 
static int close_encounter(const struct reb_simulation* const r, const struct reb_particle pi, const struct reb_particle pj)
{
    const struct reb_particle primary = r->particles[0];
    const struct reb_orbit oi = reb_orbit_from_particle(r->G, pi, primary);
    const struct reb_orbit oj = reb_orbit_from_particle(r->G, pj, primary);
    const double rhill = MAX(fabs(oi.rhill), fabs(oj.rhill));

    const double dx = pi.x - pj.x;
    const double dy = pi.y - pj.y;
    const double dz = pi.z - pj.z;
    const double d2 = dx*dx + dy*dy + dz*dz;

    return d2 <= 9.*rhill*rhill;
}

static struct reb_particle barycenter_particle(const struct hj_node *const children, const int N_children)
{
    // TODO: for now, this does not handle the case that all particles are zero-mass.
    struct reb_particle barycenter = {0};

    for (int i = 0; i < N_children; i++)
    {
        const struct reb_particle pi = children[i].barycenter_particle;
        const double m = pi.m;

        barycenter.x += m * pi.x;
        barycenter.y += m * pi.y;
        barycenter.z += m * pi.z;
        barycenter.vx += m * pi.vx;
        barycenter.vy += m * pi.vy;
        barycenter.vz += m * pi.vz;
        barycenter.ax += m * pi.ax;
        barycenter.ay += m * pi.ay;
        barycenter.az += m * pi.az;
        barycenter.m += m;
    }

    if (barycenter.m > 0.)
    {
        const double minv = 1. / barycenter.m;
        barycenter.x *= minv;
        barycenter.y *= minv;
        barycenter.z *= minv;
        barycenter.vx *= minv;
        barycenter.vy *= minv;
        barycenter.vz *= minv;
        barycenter.ax *= minv;
        barycenter.ay *= minv;
        barycenter.az *= minv;
    }

    return barycenter;
}






//----------------------------------------------------------------------------
// build tree 
static struct hj_node* reb_integrator_whfast_hj_node_add_child(struct hj_node* const node)
{
    if (node->N_children == node->N_allocated){
        const int N_allocated_new = node->N_allocated ? 2*node->N_allocated : 1;
        struct hj_node* const children_new = realloc(node->children, sizeof(struct hj_node)*N_allocated_new);
        if (children_new == NULL){
            return NULL;
        }
        node->children = children_new;
        node->N_allocated = N_allocated_new;
    }

    struct hj_node* const child = &node->children[node->N_children++];
    *child = (struct hj_node){0};
    child->particle_index = -1;
    return child;
}

static void reb_integrator_whfast_hj_node_set_leaf(struct hj_node* const node, const struct reb_particle particle, const int particle_index)
{
    reb_integrator_whfast_hj_node_free(node);
    node->barycenter_particle = particle;
    node->particle_index = particle_index;
}


static int reb_integrator_whfast_hj_build_tree_dfs(const struct reb_simulation* const r, const size_t i, struct hj_node* const cluster_node, int* const visited)
{
    struct hj_node* const leaf_node = reb_integrator_whfast_hj_node_add_child(cluster_node);
    if (leaf_node == NULL){
        return 1;
    }
    reb_integrator_whfast_hj_node_set_leaf(leaf_node, r->particles[i], (int)i);
    visited[i] = 1;

    for (size_t j=1; j<r->N; j++){
        if (!visited[j] && close_encounter(r, r->particles[i], r->particles[j])){
            int has_close_child = 0;
            for (size_t k=1; k<r->N; k++){
                if (!visited[k] && k != j && close_encounter(r, r->particles[j], r->particles[k])){
                    has_close_child = 1;
                    break;
                }
            }

            if (has_close_child){
                struct hj_node* const child_cluster_node = reb_integrator_whfast_hj_node_add_child(cluster_node);
                if (child_cluster_node == NULL){
                    return 1;
                }
                if (reb_integrator_whfast_hj_build_tree_dfs(r, j, child_cluster_node, visited)){
                    return 1;
                }
            }else{
                struct hj_node* const child_leaf_node = reb_integrator_whfast_hj_node_add_child(cluster_node);
                if (child_leaf_node == NULL){
                    return 1;
                }
                reb_integrator_whfast_hj_node_set_leaf(child_leaf_node, r->particles[j], (int)j);
                visited[j] = 1;
            }
        }
    }

    cluster_node->barycenter_particle = barycenter_particle(cluster_node->children, cluster_node->N_children);
    return 0;
}

static int reb_integrator_whfast_hj_build_tree(struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast)
{
    reb_integrator_whfast_hj_node_free(&whfast->root);

    if (r->N == 0){
        return 0;
    }

    int* const visited = calloc(r->N, sizeof(int));
    if (visited == NULL){
        return 1;
    }

    struct hj_node* const central_node = reb_integrator_whfast_hj_node_add_child(&whfast->root);
    if (central_node == NULL){
        free(visited);
        return 1;
    }
    reb_integrator_whfast_hj_node_set_leaf(central_node, r->particles[0], 0);
    visited[0] = 1;

    for (size_t i=1; i<r->N; i++){
        if (!visited[i]){
            int has_close_child = 0;
            for (size_t j=1; j<r->N; j++){
                if (!visited[j] && j != i && close_encounter(r, r->particles[i], r->particles[j])){
                    has_close_child = 1;
                    break;
                }
            }

            if (!has_close_child){
                struct hj_node* const leaf_node = reb_integrator_whfast_hj_node_add_child(&whfast->root);
                if (leaf_node == NULL){
                    free(visited);
                    reb_integrator_whfast_hj_node_free(&whfast->root);
                    return 1;
                }
                reb_integrator_whfast_hj_node_set_leaf(leaf_node, r->particles[i], (int)i);
                visited[i] = 1;
                continue;
            }

            struct hj_node* const cluster_node = reb_integrator_whfast_hj_node_add_child(&whfast->root); // check 
            if (cluster_node == NULL){
                free(visited);
                reb_integrator_whfast_hj_node_free(&whfast->root);
                return 1;
            }
            if (reb_integrator_whfast_hj_build_tree_dfs(r, i, cluster_node, visited)){
                free(visited);
                reb_integrator_whfast_hj_node_free(&whfast->root);
                return 1;
            }
        }
    }

    whfast->root.barycenter_particle = barycenter_particle(whfast->root.children, whfast->root.N_children);
    free(visited);
    return 0;
}
//----------------------------------------------------------------------------







//----------------------------------------------------------------------------
// need to be changed  
void reb_integrator_whfast_hj_from_inertial(struct reb_simulation *const r, struct reb_particle *p_jh)
{
    struct reb_particle *restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active == SIZE_MAX || r->testparticle_type == 1) ? N : r->N_active;
    



    reb_transformations_inertial_to_jacobi_posvel(particles, p_jh, particles, N, N_active);
}

void reb_integrator_whfast_hj_to_inertial(struct reb_simulation *const r, struct reb_particle *p_jh)
{
    struct reb_particle *restrict const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = (r->N_active == SIZE_MAX || r->testparticle_type == 1) ? N : r->N_active;
    reb_transformations_jacobi_to_inertial_posvel(particles, p_jh, particles, N, N_active);
}

/***************************** 
 * Interaction Hamiltonian  */
void reb_integrator_whfast_hj_interaction_step(struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    const double G = r->G;
    struct reb_particle* particles = r->particles;
    const double m0 = particles[0].m;
    reb_transformations_inertial_to_jacobi_acc(particles, p_jh, particles, N, N_active);
    double eta = m0;
    for (size_t i=1;i<N;i++){
        // Eq 132
        const struct reb_particle pji = p_jh[i];
        if (i<N_active){
            eta += pji.m;
        }
        // ax was calculate by update_acceleration O(N^2), last term from the right
        p_jh[i].vx += _dt * pji.ax;
        p_jh[i].vy += _dt * pji.ay;
        p_jh[i].vz += _dt * pji.az;
        
        // Additional Jacobi terms (second term from the right)
        const double rj2i = 1./(pji.x*pji.x + pji.y*pji.y + pji.z*pji.z);
        const double rji  = sqrt(rj2i);
        const double rj3iM = rji*rj2i*G*eta;
        const double prefac1 = _dt*rj3iM;
        p_jh[i].vx += prefac1*pji.x;
        p_jh[i].vy += prefac1*pji.y;
        p_jh[i].vz += prefac1*pji.z;
    }
}

/***************************** 
 * DKD Scheme                */

void reb_integrator_whfast_hj_kepler_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    const double m0 = r->particles[0].m;
    const double G = r->G;
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type ==1)?N:r->N_active;
    double eta = m0;
    for (size_t i=1;i<N;i++){
        if (i<N_active){
            eta += p_jh[i].m;
        }
        reb_integrator_whfast_kepler_solver(&p_jh[i], eta*G, _dt, r);
    }
}

void reb_integrator_whfast_hj_com_step(const struct reb_simulation* const r, struct reb_particle* p_jh, const double _dt){
    (void)r;
    p_jh[0].x += _dt*p_jh[0].vx;
    p_jh[0].y += _dt*p_jh[0].vy;
    p_jh[0].z += _dt*p_jh[0].vz;
}
//----------------------------------------------------------------------------



void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_hj_state* whfast = state;
    if (reb_integrator_whfast_hj_build_tree(r, whfast)){
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for the tree.");
        return;
    }

    // TODO: transform inertial coordinates into the HJ tree.
    // TODO: apply the recursive H_A half step.
    // TODO: transform tree coordinates back to inertial coordinates.
    // TODO: update accelerations and apply the recursive H_B full step.
    // TODO: apply the recursive H_A half step.
    // TODO: transform tree coordinates back to inertial coordinates.
    // TODO: advance r->t and r->dt_last_done.

    /*
     * Reference: old flat Jacobi WHFast-HJ step.
     *
     * const double dt = r->dt;
     *
     * // H_A half step
     * reb_integrator_whfast_hj_from_inertial(r, whfast->p_jh);
     * reb_integrator_whfast_hj_kepler_step(r, whfast->p_jh, r->dt/2.);
     * reb_integrator_whfast_hj_com_step(r, whfast->p_jh, r->dt/2.);
     * reb_integrator_whfast_hj_to_inertial(r, whfast->p_jh);
     *
     * // H_B full step
     * reb_simulation_update_acceleration(r);
     * reb_integrator_whfast_hj_interaction_step(r, whfast->p_jh, dt);
     *
     * // H_A half step
     * reb_integrator_whfast_hj_kepler_step(r, whfast->p_jh, r->dt/2.);
     * reb_integrator_whfast_hj_com_step(r, whfast->p_jh, r->dt/2.);
     *
     * reb_integrator_whfast_hj_to_inertial(r, whfast->p_jh);
     *
     * //   0 1 2 3 4 5
     * // 0   x x x x x
     * // 1 x   x x x x
     * // 2 x x   x x x
     * // 3 x x x   x x
     * // 4 x x x x   x
     * // 5 x x x x x
     * r->t += r->dt;
     * r->dt_last_done = r->dt;
     */

    reb_simulation_error(r, "WHFast HJ tree stepping is not implemented yet.");
    r->status = REB_STATUS_GENERIC_ERROR;
}
