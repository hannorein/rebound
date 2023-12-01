# Naming Convention

Starting with REBOUND version 4.0 we try to keep all variable, structure, class, and function names adhere to a naming convention which is outlined in this document.
We do this because this will make it easier to users to understand what a variable or function does.
There is some tension between the C and python side but this document provides a clear map to translate a name from c to python and vice versa.


- On the python side, we should follow [PEP8](https://peps.python.org/pep-0008/#function-and-variable-names). 
- For both python and C, function names and variables use lowercase, with words separated by underscores as necessary to improve readability.
- Structure names in C start with the prefix `reb_` followed by small letters, words are separated by underscores as necessary to improve readability, e.g. `reb_hash_pointer_pair `.
- Structures in c do not use typedef. They are always referred to `struct`, e.g. `struct reb_simulation`.
- Python class names use the `CapWords` convention and ignore the `reb_` prefix of the corresponding c structure. For example: `reb_simulation` <-> `Simulation`, `reb_simulationarchive` <-> `Simulationarchive`, `reb_hash_pointer_pair` <-> `HashPointerPair`
- Function names in C start with the main object they operate on. For example `reb_simulation_move_to_com()` operates on the `reb_simulation` object.
- These functions often translate to instance method in python. In that case the object the function operates on is class object. For example the C function `reb_simulation_move_to_com()` becomes the instance method `rebound.Simulation.move_to_com()` in python. 
- Function names only include a verb if the function changes the state of an object. For example: `reb_simulation_move_to_com` changes the simulation by moving the center of mass frame and therefore includes the verb move. `reb_simulation_com` on the other hand simply calculates the center of mass and does not change the simulation. The name does therefore not include a verb. There are cases where it does not make sense to enforce this rule. 
- Function names try to avoid using "set" and "get". For example, instead of `reb_get_com()` we use `reb_simulation_com()`.
- Variables describing memory allocation counts have names starting with `N_allocated`. For example: `N_allocated_collisions`.
- Functions that are related to memory management use the following verbs: `create` allocates and initializes an object. `free` frees the memory of an object (and all the objects it owns). `init` does not allocated an object itself, but merely initializes it with default values. 
- As with all functions those related to memory allocation also start with the object they operate on, then followed by the verb. For example: `reb_simulation_create()`.
- Simulationarchive is one word. All functions related to a Simulationarchive in C use `_simulationarchive_` and `Simulationarchive` in python (not `SimulationArchive`). 
