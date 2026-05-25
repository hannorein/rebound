
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
    if (node == NULL){
        return;
    }
    reb_integrator_whfast_hj_node_free(node->primary);
    reb_integrator_whfast_hj_node_free(node->secondary);
    free(node);
}

void reb_integrator_whfast_hj_free(void* p){
    struct reb_integrator_whfast_hj_state* whfast = p;
    reb_integrator_whfast_hj_node_free(whfast->root);
    free(whfast);
}

static struct hj_node* reb_integrator_whfast_hj_node_create_leaf(const struct reb_particle particle, const int particle_index)
{
    struct hj_node* const node = calloc(1, sizeof(struct hj_node));
    if (node == NULL){
        return NULL;
    }
    node->barycenter_particle = particle;
    node->particle_index = particle_index;
    return node;
}

static struct hj_node* reb_integrator_whfast_hj_node_create_binary(struct hj_node* const a, struct hj_node* const b)
{
    struct hj_node* const node = calloc(1, sizeof(struct hj_node));
    if (node == NULL){
        return NULL;
    }

    if (a->barycenter_particle.m >= b->barycenter_particle.m){
        node->primary = a;
        node->secondary = b;
    }else{
        node->primary = b;
        node->secondary = a;
    }
    node->barycenter_particle = reb_particle_com_of_pair(node->primary->barycenter_particle, node->secondary->barycenter_particle);
    node->particle_index = -1;
    return node;
}

static void reb_integrator_whfast_hj_node_list_free(struct hj_node** const nodes, const size_t N)
{
    for (size_t i=0; i<N; i++){
        reb_integrator_whfast_hj_node_free(nodes[i]);
    }
    free(nodes);
}

static struct reb_orbit reb_integrator_whfast_hj_orbit_for_pair(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    if (a->barycenter_particle.m >= b->barycenter_particle.m){
        return reb_orbit_from_particle(G, b->barycenter_particle, a->barycenter_particle);
    }
    return reb_orbit_from_particle(G, a->barycenter_particle, b->barycenter_particle);
}

static double reb_integrator_whfast_hj_bound_pair_timescale(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    const struct reb_orbit o = reb_integrator_whfast_hj_orbit_for_pair(G, a, b);
    if (isfinite(o.P) && o.P > 0.){
        return o.P;
    }
    return INFINITY;
}

static double reb_integrator_whfast_hj_unbound_pair_timescale(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    const struct reb_orbit o = reb_integrator_whfast_hj_orbit_for_pair(G, a, b);
    const double mu = G*(a->barycenter_particle.m + b->barycenter_particle.m);
    const double denominator = mu*(1.0 + o.e);
    if (!isfinite(o.e) || denominator <= 0.){
        return INFINITY;
    }

    const double q = fabs(o.e*(1.0 - o.e));
    const double Pperi = sqrt(q*q*q/denominator);
    if (isfinite(Pperi)){
        return Pperi;
    }
    return INFINITY;
}

//----------------------------------------------------------------------------
// build binary orbit hierarchy tree
static int reb_integrator_whfast_hj_build_tree(struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast)
{
    reb_integrator_whfast_hj_node_free(whfast->root);
    whfast->root = NULL;

    size_t N = r->N;
    if (N == 0){
        return 0;
    }

    struct hj_node** const nodes = calloc(N, sizeof(struct hj_node*));
    if (nodes == NULL){
        return 1;
    }

    for (size_t i=0; i<N; i++){
        nodes[i] = reb_integrator_whfast_hj_node_create_leaf(r->particles[i], (int)i);
        if (nodes[i] == NULL){
            reb_integrator_whfast_hj_node_list_free(nodes, i);
            return 1;
        }
    }

    while (N > 1){
        size_t imin = SIZE_MAX;
        size_t jmin = SIZE_MAX;
        double Pmin = INFINITY;

        for (size_t i=1; i<N; i++){
            for (size_t j=0; j<i; j++){
                const double P = reb_integrator_whfast_hj_bound_pair_timescale(r->G, nodes[i], nodes[j]);
                if (P < Pmin){
                    Pmin = P;
                    imin = i;
                    jmin = j;
                }
            }
        }

        if (Pmin == INFINITY){
            for (size_t i=1; i<N; i++){
                for (size_t j=0; j<i; j++){
                    const double Pperi = reb_integrator_whfast_hj_unbound_pair_timescale(r->G, nodes[i], nodes[j]);
                    if (Pperi < Pmin){
                        Pmin = Pperi;
                        imin = i;
                        jmin = j;
                    }
                }
            }
        }

        if (Pmin == INFINITY){
            reb_integrator_whfast_hj_node_list_free(nodes, N);
            return 1;
        }

        struct hj_node* const merged = reb_integrator_whfast_hj_node_create_binary(nodes[imin], nodes[jmin]);
        if (merged == NULL){
            reb_integrator_whfast_hj_node_list_free(nodes, N);
            return 1;
        }

        nodes[imin] = merged;
        nodes[jmin] = nodes[N-1];
        N--;
    }

    whfast->root = nodes[0];
    free(nodes);
    return 0;
}

//----------------------------------------------------------------------------
// Tree-Jacobi coordinate transforms.
static int reb_integrator_whfast_hj_node_is_leaf(const struct hj_node* const node)
{
    return node != NULL && node->primary == NULL && node->secondary == NULL;
}

void reb_integrator_whfast_hj_from_inertial(struct reb_simulation* const r, struct hj_node* const node)
{
    if (node == NULL){
        return;
    }

    if (reb_integrator_whfast_hj_node_is_leaf(node)){
        if (node->particle_index >= 0 && (size_t)node->particle_index < r->N){
            node->barycenter_particle = r->particles[node->particle_index];
        }
        node->jacobi_particle = (struct reb_particle){0};
        return;
    }

    reb_integrator_whfast_hj_from_inertial(r, node->primary);
    reb_integrator_whfast_hj_from_inertial(r, node->secondary);

    const struct reb_particle primary = node->primary->barycenter_particle;
    const struct reb_particle secondary = node->secondary->barycenter_particle;

    if (primary.m + secondary.m > 0.){
        node->barycenter_particle = reb_particle_com_of_pair(primary, secondary);
    }else{
        node->barycenter_particle = primary;
        node->barycenter_particle.m = 0.;
    }

    node->jacobi_particle = (struct reb_particle){0};
    node->jacobi_particle.m = secondary.m;
    node->jacobi_particle.x = secondary.x - primary.x;
    node->jacobi_particle.y = secondary.y - primary.y;
    node->jacobi_particle.z = secondary.z - primary.z;
    node->jacobi_particle.vx = secondary.vx - primary.vx;
    node->jacobi_particle.vy = secondary.vy - primary.vy;
    node->jacobi_particle.vz = secondary.vz - primary.vz;
    node->jacobi_particle.ax = secondary.ax - primary.ax;
    node->jacobi_particle.ay = secondary.ay - primary.ay;
    node->jacobi_particle.az = secondary.az - primary.az;
}

static void reb_integrator_whfast_hj_reconstruct_child_barycenters(struct hj_node* const node)
{
    struct reb_particle primary = node->primary->barycenter_particle;
    struct reb_particle secondary = node->secondary->barycenter_particle;
    const struct reb_particle barycenter = node->barycenter_particle;
    const struct reb_particle jacobi = node->jacobi_particle;

    const double M_primary = primary.m;
    const double M_secondary = secondary.m;
    const double M_total = M_primary + M_secondary;

    if (M_total > 0.){
        const double primary_offset = M_secondary/M_total;
        const double secondary_offset = M_primary/M_total;

        primary.x = barycenter.x - primary_offset*jacobi.x;
        primary.y = barycenter.y - primary_offset*jacobi.y;
        primary.z = barycenter.z - primary_offset*jacobi.z;
        primary.vx = barycenter.vx - primary_offset*jacobi.vx;
        primary.vy = barycenter.vy - primary_offset*jacobi.vy;
        primary.vz = barycenter.vz - primary_offset*jacobi.vz;
        primary.ax = barycenter.ax - primary_offset*jacobi.ax;
        primary.ay = barycenter.ay - primary_offset*jacobi.ay;
        primary.az = barycenter.az - primary_offset*jacobi.az;

        secondary.x = barycenter.x + secondary_offset*jacobi.x;
        secondary.y = barycenter.y + secondary_offset*jacobi.y;
        secondary.z = barycenter.z + secondary_offset*jacobi.z;
        secondary.vx = barycenter.vx + secondary_offset*jacobi.vx;
        secondary.vy = barycenter.vy + secondary_offset*jacobi.vy;
        secondary.vz = barycenter.vz + secondary_offset*jacobi.vz;
        secondary.ax = barycenter.ax + secondary_offset*jacobi.ax;
        secondary.ay = barycenter.ay + secondary_offset*jacobi.ay;
        secondary.az = barycenter.az + secondary_offset*jacobi.az;
    }else{
        primary.x = barycenter.x;
        primary.y = barycenter.y;
        primary.z = barycenter.z;
        primary.vx = barycenter.vx;
        primary.vy = barycenter.vy;
        primary.vz = barycenter.vz;
        primary.ax = barycenter.ax;
        primary.ay = barycenter.ay;
        primary.az = barycenter.az;

        secondary.x = barycenter.x + jacobi.x;
        secondary.y = barycenter.y + jacobi.y;
        secondary.z = barycenter.z + jacobi.z;
        secondary.vx = barycenter.vx + jacobi.vx;
        secondary.vy = barycenter.vy + jacobi.vy;
        secondary.vz = barycenter.vz + jacobi.vz;
        secondary.ax = barycenter.ax + jacobi.ax;
        secondary.ay = barycenter.ay + jacobi.ay;
        secondary.az = barycenter.az + jacobi.az;
    }

    node->primary->barycenter_particle = primary;
    node->secondary->barycenter_particle = secondary;
}

void reb_integrator_whfast_hj_to_inertial(struct reb_simulation* const r, struct hj_node* const node)
{
    if (node == NULL){
        return;
    }

    if (reb_integrator_whfast_hj_node_is_leaf(node)){
        if (node->particle_index >= 0 && (size_t)node->particle_index < r->N){
            struct reb_particle* const particle = &r->particles[node->particle_index];
            const struct reb_particle source = node->barycenter_particle;
            particle->x = source.x;
            particle->y = source.y;
            particle->z = source.z;
            particle->vx = source.vx;
            particle->vy = source.vy;
            particle->vz = source.vz;
        }
        return;
    }

    reb_integrator_whfast_hj_reconstruct_child_barycenters(node);
    reb_integrator_whfast_hj_to_inertial(r, node->primary);
    reb_integrator_whfast_hj_to_inertial(r, node->secondary);
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

    reb_simulation_error(r, "WHFast HJ tree stepping is not implemented yet.");
    r->status = REB_STATUS_GENERIC_ERROR;
}
