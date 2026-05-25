/**
 * Print-only diagnostic for a three-particle tree-Jacobi hierarchy.
 *
 * The particles are arranged so the tree builder should form:
 *
 *     (1, (2, 3))
 *
 * This program intentionally does not assert the tree shape. It prints the
 * inertial particles, the tree structure and node masses, the HJ coordinates
 * after from_inertial(), and the inertial particles after to_inertial().
 */

#include <math.h>
#include <stdio.h>

#include "rebound.h"
#include "../../../../src/integrator_whfast_HJ.c"

#define N_PARTICLES 3

static int particle_label(const int particle_index)
{
    return particle_index + 1;
}

static void print_indent(const int depth)
{
    for (int i = 0; i < depth; i++)
    {
        printf("  ");
    }
}

static void print_particle_fields(const struct reb_particle p)
{
    printf("m=% .5g  r=(% .5g, % .5g, % .5g)  v=(% .5g, % .5g, % .5g)",
           p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz);
}

static void print_particle_line(const int label, const struct reb_particle p)
{
    printf("particle %d: ", label);
    print_particle_fields(p);
    printf("\n");
}

static void print_particles(const struct reb_simulation *const r, const char *const title)
{
    printf("\n%s\n", title);
    for (int i = 0; i < r->N; i++)
    {
        print_particle_line(particle_label(i), r->particles[i]);
    }
}

static int is_leaf_node(const struct hj_node *const node)
{
    return node != NULL && node->primary == NULL && node->secondary == NULL;
}

static void print_tree_structure_node(const struct hj_node *const node, const int depth, const char *const branch)
{
    if (node == NULL)
    {
        print_indent(depth);
        printf("%s: NULL\n", branch);
        return;
    }

    print_indent(depth);
    if (is_leaf_node(node))
    {
        printf("%s leaf particle %d: ", branch, particle_label(node->particle_index));
        print_particle_fields(node->barycenter_particle);
        printf("\n");
        return;
    }

    printf("%s internal node: ", branch);
    print_particle_fields(node->barycenter_particle);
    printf("\n");
    print_tree_structure_node(node->primary, depth + 1, "primary");
    print_tree_structure_node(node->secondary, depth + 1, "secondary");
}

static void print_tree_structure(const struct hj_node *const root)
{
    printf("\nTree structure and mass at each node\n");
    print_tree_structure_node(root, 0, "root");
}

static void print_hj_node(const struct hj_node *const node, const int depth, const char *const branch)
{
    if (node == NULL)
    {
        print_indent(depth);
        printf("%s: NULL\n", branch);
        return;
    }

    print_indent(depth);
    if (is_leaf_node(node))
    {
        printf("%s leaf particle %d barycenter: ", branch, particle_label(node->particle_index));
        print_particle_fields(node->barycenter_particle);
        printf("\n");
        return;
    }

    printf("%s internal node barycenter: ", branch);
    print_particle_fields(node->barycenter_particle);
    printf("\n");

    print_indent(depth);
    printf("%s internal node HJ coordinate: ", branch);
    print_particle_fields(node->jacobi_particle);
    printf("\n");

    print_hj_node(node->primary, depth + 1, "primary");
    print_hj_node(node->secondary, depth + 1, "secondary");
}

static void print_hj_coordinates(const struct hj_node *const root)
{
    printf("\nHJ coordinates after from_inertial()\n");
    print_hj_node(root, 0, "root");
}

static void add_particle(struct reb_simulation *const r, const struct reb_particle p)
{
    reb_simulation_add(r, p);
}

static struct reb_simulation *setup_three_particle_tree_simulation(void)
{
    struct reb_simulation *const r = reb_simulation_create();
    r->G = 1.0;

    const double m1 = 1.0;
    const double m2 = 0.02;
    const double m3 = 0.01;
    const double m23 = m2 + m3;
    const double pair_separation = 0.02;
    const double pair_center_x = 5.0;
    const double pair_relative_speed = sqrt(r->G * m23 / pair_separation);

    struct reb_particle p1 = {0};
    p1.m = m1;
    p1.x = 0.0;
    p1.y = 0.0;
    p1.z = 0.0;

    struct reb_particle p2 = {0};
    p2.m = m2;
    p2.x = pair_center_x - (m3 / m23) * pair_separation;
    p2.y = 0.0;
    p2.z = 0.0;
    p2.vx = 0.0;
    p2.vy = -(m3 / m23) * pair_relative_speed;
    p2.vz = 0.0;

    struct reb_particle p3 = {0};
    p3.m = m3;
    p3.x = pair_center_x + (m2 / m23) * pair_separation;
    p3.y = 0.0;
    p3.z = 0.0;
    p3.vx = 0.0;
    p3.vy = (m2 / m23) * pair_relative_speed;
    p3.vz = 0.0;

    add_particle(r, p1);
    add_particle(r, p2);
    add_particle(r, p3);

    return r;
}

int main(void)
{
    printf("three particle (1,(2,3))\n");

    struct reb_simulation *const r = setup_three_particle_tree_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    print_particles(r, "Initial inertial coordinates and masses");

    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0)
    {
        fprintf(stderr, "Failed to build HJ tree.\n");
        reb_simulation_free(r);
        return 1;
    }

    print_tree_structure(whfast.root);

    reb_integrator_whfast_hj_from_inertial(r, whfast.root);
    print_hj_coordinates(whfast.root);

    reb_integrator_whfast_hj_to_inertial(r, whfast.root);
    print_particles(r, "Inertial coordinates after to_inertial()");

    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return 0;
}
