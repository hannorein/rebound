# Collisions

## Detecting collisions

REBOUND comes with several collision detection modules. 
These modules check for physical collisions (the distance between two particles is closer than the sum of the radii), not close encounters.
For a collision to occur between two particles, at least one of them needs to have a finite radius and collision detection needs to be turned on (it is turned off by default).


### No collisions
By default REBOUND does not search for collisions. 
You can manually set the collision routine to NONE with the following code:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_NONE;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.collision = "none"
    ```

### Direct
The direct collision detection module is a brute force collision search and scales as $O(N^2)$.
It checks for instantaneous overlaps between every particle pair. 
The following code enables this module:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_DIRECT;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.collision = "direct"
    ```

!!! Important
    This method checks for instantaneous overlaps. It does this only after each timestep.
    This means that if the timestep is large enough for particles to pass completely through each other, then the collision will be missed. 
    


### Line
This is a brute force collision search and scales as $O(N^2)$ but compared to the direct method described above, this algorithm checks for overlapping particles during the timestep (not just at the end).
It assumes particles travelled along straight lines during the timestep and might therefore miss some collisions.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_LINE;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.collision = "line"
    ```

### Tree
This method uses an oct-tree to check for overlapping particles at the end of the timestep.
When a large number of particles $N$ is used, this method scales as $O(N log(N))$, rather than $O(N^2)$ for the direct search.
Note that you need to initialize the simulation box whenever you want to use the tree.
Below is an example on how to enable the tree based collision search.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_configure_box(r, 10, 1, 1, 1); # confine the simulation to a box of size 10
    r->collision = REB_COLLISION_TREE;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.configure_box(10)   # confine the simulation to a box of size 10
    sim.collision = "tree"
    ```


### Linetree
Similar to the tree method, this method also uses an oct-tree and has a scaling of $O(N log(N))$.  
It checks for overlapping trajectories during the last timestep, not only for overlapping particles at the end of the timestep.
It might still miss some collisions because it assumes that particles travel along straight lines.


Below is an example on how to enable the line-tree collision search.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_configure_box(r, 10, 1, 1, 1); # confine the simulation to a box of size 10
    r->collision = REB_COLLISION_LINETREE;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.configure_box(10)   # confine the simulation to a box of size 10
    sim.collision = "linetree"
    ```

## Resolving collisions

Once a collision has been detected, you have a choice on what to do next.
You might just want to merge particles, let them bounce off each other, or simply keep a log of all collisions that occurred. 

REBOUND comes with several built-in collision resolve functions. 
You can also write your own.

Internally this functionality is implemented using a [function pointer](https://www.cprogramming.com/tutorial/function-pointers.html). 
You can set this pointer to a function that should be called when a collision occurs, whether it be a built-in function or your own. 

### Halt

This function resolves a collision by simply halting the integration and setting the `status` flag in the simulation to `REB_STATUS_COLLISION`. 
In python this will raise the `Collision` exception.
This is the default.
It can also be set manually using the following syntax:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_halt;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.collision = "direct"
    sim.collision_resolve = "halt"
    ```

### Hard-sphere

This assumes a hard-sphere collision and uses the `coefficient_of_restitution` function pointer in `struct reb_simulation` to determine coefficient of restitution which can be velocity dependent
It conserves momentum and mass.
Depending on the coefficient of restitution, it also conserves energy.

The following example shows how to set up a hard-sphere collision resolve function and a direct collision detection routine.

=== "C"
    ```c
    double coefficient_of_restitution_constant(const struct reb_simulation* const r, double v){
        // v is the normal impact velocity.
        // Here, we just use a constant coefficient of restitution
        return 0.5;
    }
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_DIRECT;
    r->coefficient_of_restitution = coefficient_of_restitution_constant;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    ```

=== "Python"
    ```python
    def coefficient_of_restitution_constant(r, v):
        # v is the normal impact velocity.
        # Here, we just use a constant coefficient of restitution
        return 0.5
    sim = rebound.Simulation()
    sim.collision = "direct"
    sim.coefficient_of_restitution = coefficient_of_restitution_constant
    sim.collision_resolve = "hardsphere"
    ```


### Merge
 
This function merges the two colliding particles. 
It conserves mass, momentum and volume, but not energy.  
The particle with the higher index will be removed. 

The following example shows how to set up a hard-sphere collision resolve function and a direct collision detection routine.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.collision = "direct"
    sim.collision_resolve = "merge"
    ```

### Custom function
You can write your own collision resolve function. 
In your function, you can update the properties of the particles involved in the collision.
The return value of your function determines if a particle gets removed.

- `0`: don't remove either particle from the simulation
- `1`: remove the first particle (`p1`) from the simulation
- `2`: remove the second particle (`p2`) from the simulation
- `3`: remove both particles from the simulation

Here is a short example on how to write a simple custom collision resolve function:

=== "C"
    ```c
    int collision_print_only(struct reb_simulation* const r, struct reb_collision c){
        printf("%f\t", r->p);
        printf("%f\t", r->particles[c.p1].x);    // x position of particle 1
        printf("%f\n", r->particles[c.p2].x);    // x position of particle 2
        return 0; // Don't remove either particle
    }

    int main(int argc, char* argv[]){
        struct reb_simulation* r = reb_simulation_create();
        r->collision = REB_COLLISION_DIRECT;
        r->collision_resolve = collision_print_only;
    }
    ```

=== "Python"
    ```python
    def collision_print_only(sim_pointer, collision):
        sim = sim_pointer.contents           # get simulation object from pointer
        print(sim.t)                         # print time 
        print(sim.particles[collision.p1].x) # x position of particle 1
        print(sim.particles[collision.p2].x) # x position of particle 2
        return 0                             # Don't remove either particle
    
    sim = rebound.Simulation()
    sim.collision = "direct"
    sim.collision_resolve = collision_print_only
    ```
The first argument of the collision resolve function is a pointer to the simulation. 
The second argument is a `reb_collision` structure. 
It contains information about which particles are involved in the collision and, for periodic or shear-periodic [boundary conditions](boundaryconditions.md), if the collision occurred across a boundary:

`int p1`
:   Index corresponding to one of the colliding particles

`int p2`
:   Index corresponding to one of the colliding particles

`struct reb_vec6d gb`
:   Shift of particle p1 due to a collision across periodic and shearing sheet boundaries. All entries are zero if a normal collision occurs.
