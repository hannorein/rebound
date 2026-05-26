
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

#define WHFAST_HJ_UNBOUND_DISTANCE_FACTOR 0.3

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
    const double M = a->barycenter_particle.m + b->barycenter_particle.m;
    const double mu = G*M;
    if (!isfinite(mu) || mu <= 0. || !isfinite(o.e) || o.e <= 1.){
        return INFINITY;
    }

    const struct reb_particle pa = a->barycenter_particle;
    const struct reb_particle pb = b->barycenter_particle;

    const double rx = pb.x - pa.x;
    const double ry = pb.y - pa.y;
    const double rz = pb.z - pa.z;
    const double vx = pb.vx - pa.vx;
    const double vy = pb.vy - pa.vy;
    const double vz = pb.vz - pa.vz;

    const double hx = ry*vz - rz*vy;
    const double hy = rz*vx - rx*vz;
    const double hz = rx*vy - ry*vx;
    const double h2 = hx*hx + hy*hy + hz*hz;
    const double denominator = mu*(1.0 + o.e);
    if (!isfinite(h2) || h2 <= 0. || !isfinite(denominator) || denominator <= 0.){
        return INFINITY;
    }

    const double rp = h2/denominator;
    if (!isfinite(rp) || rp <= 0.){
        return INFINITY;
    }

    const double tperi = sqrt(rp*rp*rp/denominator);
    return (isfinite(tperi) && tperi > 0.) ? tperi : INFINITY;
}

static double reb_integrator_whfast_hj_pair_distance(const struct hj_node* const a, const struct hj_node* const b)
{
    const double dx = b->barycenter_particle.x - a->barycenter_particle.x;
    const double dy = b->barycenter_particle.y - a->barycenter_particle.y;
    const double dz = b->barycenter_particle.z - a->barycenter_particle.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
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
        double dmin_bound = INFINITY;

        for (size_t i=1; i<N; i++){
            for (size_t j=0; j<i; j++){
                const double P = reb_integrator_whfast_hj_bound_pair_timescale(r->G, nodes[i], nodes[j]);
                if (P < Pmin){
                    Pmin = P;
                    imin = i;
                    jmin = j;
                    dmin_bound = reb_integrator_whfast_hj_pair_distance(nodes[i], nodes[j]);
                }
            }
        }

        for (size_t i=1; i<N; i++){
            for (size_t j=0; j<i; j++){
                const double Pperi = reb_integrator_whfast_hj_unbound_pair_timescale(r->G, nodes[i], nodes[j]);
                if (Pperi < Pmin && reb_integrator_whfast_hj_pair_distance(nodes[i], nodes[j]) < WHFAST_HJ_UNBOUND_DISTANCE_FACTOR*dmin_bound){
                    Pmin = Pperi;
                    imin = i;
                    jmin = j;
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
    const double M_primary = primary.m;
    const double M_secondary = secondary.m;
    const double M_total = M_primary + M_secondary;

    if (M_total > 0.){
        node->barycenter_particle = reb_particle_com_of_pair(primary, secondary);
    }else{
        node->barycenter_particle = primary;
        node->barycenter_particle.m = 0.;
    }

    node->jacobi_particle = (struct reb_particle){0};
    node->jacobi_particle.m = (M_total > 0.) ? M_primary*M_secondary/M_total : 0.;
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
static void reb_integrator_whfast_hj_interaction_step_node(const struct reb_simulation* const r, struct hj_node* const node, const double _dt){
    if (node == NULL || reb_integrator_whfast_hj_node_is_leaf(node)){
        return;
    }

    struct reb_particle* const p = &node->jacobi_particle;
    const double eta = node->primary->barycenter_particle.m + node->secondary->barycenter_particle.m;
    p->vx += _dt*p->ax;
    p->vy += _dt*p->ay;
    p->vz += _dt*p->az;

    const double rj2i = 1./(p->x*p->x + p->y*p->y + p->z*p->z);
    const double rji = sqrt(rj2i);
    const double prefac = _dt*r->G*eta*rji*rj2i;
    p->vx += prefac*p->x;
    p->vy += prefac*p->y;
    p->vz += prefac*p->z;

    reb_integrator_whfast_hj_interaction_step_node(r, node->primary, _dt);
    reb_integrator_whfast_hj_interaction_step_node(r, node->secondary, _dt);
}

void reb_integrator_whfast_hj_interaction_step(struct reb_simulation* const r, struct hj_node* const root, const double _dt){
    reb_integrator_whfast_hj_from_inertial(r, root);
    reb_integrator_whfast_hj_interaction_step_node(r, root, _dt);
}

/***************************** 
 * DKD Scheme                */

void reb_integrator_whfast_hj_kepler_step(const struct reb_simulation* const r, struct hj_node* const node, const double _dt){
    if (node == NULL || reb_integrator_whfast_hj_node_is_leaf(node)){
        return;
    }

    const double eta = node->primary->barycenter_particle.m + node->secondary->barycenter_particle.m;
    reb_integrator_whfast_kepler_solver(&node->jacobi_particle, eta*r->G, _dt, r);
    reb_integrator_whfast_hj_kepler_step(r, node->primary, _dt);
    reb_integrator_whfast_hj_kepler_step(r, node->secondary, _dt);
}

void reb_integrator_whfast_hj_com_step(const struct reb_simulation* const r, struct hj_node* const root, const double _dt){
    (void)r;
    if (root == NULL){
        return;
    }
    root->barycenter_particle.x += _dt*root->barycenter_particle.vx;
    root->barycenter_particle.y += _dt*root->barycenter_particle.vy;
    root->barycenter_particle.z += _dt*root->barycenter_particle.vz;
}
//----------------------------------------------------------------------------

void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_hj_state* whfast = state;
    const double dt = r->dt;
    if (reb_integrator_whfast_hj_build_tree(r, whfast)){
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for the tree.");
        return;
    }

    reb_integrator_whfast_hj_from_inertial(r, whfast->root);

    reb_integrator_whfast_hj_kepler_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_com_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast->root);

    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    reb_simulation_update_acceleration(r);
    reb_integrator_whfast_hj_interaction_step(r, whfast->root, dt);

    reb_integrator_whfast_hj_kepler_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_com_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast->root);

    r->t += dt;
    r->dt_last_done = dt;
}
