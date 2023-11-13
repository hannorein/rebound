# Gravity solvers

## Basic
`REB_GRAVITY_BASIC`

The basic gravity routine works is the default. It works in most cases. 
It uses direct summation to calculate gravitational forces between all particle pairs.
OpenMP parallelization is implemented. If OpenMP is turned on, the scaling is $O(N^2)$, otherwise, it is $O(\frac12 N^2)$, where $N$ is the number of particles. 

## Compensated
`REB_GRAVITY_COMPENSATED`

This routine also uses direct summation but in addition makes use of compensated summation to minimize round-off errors. 
There are only a few special cases where the round-off error in force calculations has a dominant effect. In most cases, the basic gravity routine is faster and equally accurate.

## Tree
`REB_GRAVITY_TREE`          

This method uses an oct tree (Barnes and Hut 1986) to approximate self-gravity. It scales as  $O(N \log(N))$.

## Tree
`REB_GRAVITY_JACOBI`        

Direct summation, scales as $O(N^2)$, includes special terms needed for some symplectic integrators.

## None
`REB_GRAVITY_NONE`          

By using this gravity routine, no self-gravity calculated. It is still possible to include additional forces. 
