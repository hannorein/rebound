# Binary Format

REBOUND comes with its own binary format.
The binary format allows you to store a current simulation state to a file or to memory.
The binary format is also used when you make a copy of a simulation or when you compare two simulations with each other.
The Simulationarchive is an extension of the binary format which allows you to store multiple snapshots of a simulation in one file.
This page explains the details of the binary format.
It is mainly intended for people who wish to extend the built-in REBOUND functionality.
You do not need to know those details if you're only working with binary files to save and load simulations.

REBOUND uses two structures for the binary files:

```c
struct reb_binary_field {
    uint32_t type; 
    uint64_t size;
};
```

and 

```c
struct reb_simulationarchive_blob {
    int32_t index;
    int32_t offset_prev;
    int32_t offset_next;
};
```

!!! note
    Before version 3.18, the offset datatype was `int16_t`. This caused problems for simulations with a large number of particles and has since been change to `int32_t`.

## Binary file (one snapshot)
You create a binary file if you save a simulation
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_simulation_save_to_file(r, "snapshot.bin");
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    // ... setup simulation ...
    sim.save_to_file("snapshot.bin")
    ```
Such a binary file with one snapshot is simply a set of `reb_binaryfield`s followed by one `reb_simulationarchive_blob` at the end, for example:

```
reb_binary_field:
    type: DT
    size: 8 bytes

8 bytes of data representing the value of DT

reb_binary_field:
    type: PARTICLES
    size: 128 bytes

128 bytes of data representing the values of PARTICLES

...

reb_binary_field:
    type: END
    size: 0

reb_simulationarchive_blob:
    index: 0
    offset_prev: 0
    offset_next: 0
```

Each of the binary fields provides the context (type and size) for the data that immediately follows the field.
The type is an integer defined in the `reb_binary_field_descriptor_list` (see below).
The last binary field of type `9999` (`end`) to indicate that the snapshot ends here. 

!!! note
    Before version 3.27 data was encoded using the enum `REB_BINARY_FIELD_TYPE` instead of `reb_binary_field_descriptor_list`.


## Simulationarchive file (multiple snapshots)
The binary file above can also be interpreted as a Simulationarchive with one snapshot. 
You can append many (millions!) of snapshots to a binary file.
REBOUND only stores data that has changed since the original snapshot (typically the particle data, time, etc).
This allows for a very compact file size, while still maintaining bit-wise reproducibility. 

Each snapshot is separated by a `reb_simulationarchive_blob`. 
The blob contains the offset to the previous and next blobs. 
This allows REBOUND to quickly jump from one blob in the archive to the next.
Between the blobs are the same `reb_binary_field`s we already encountered for a binary file with one snapshot.
Thus, a Simulationarchive file with multiple snapshots looks something like this:

```
reb_binary_field:
    type: DT
    size: 8 bytes

8 bytes of data representing the value of DT

... more reb_binary_fields ...

reb_binary_field:
    type: END
    size: 0

reb_simulationarchive_blob:
    index: 0
    offset_prev: 0
    offset_next: 256 (offset to the next blob)

reb_binary_field:
    type: DT
    size: 8 bytes

8 bytes of data representing the value of DT

... more reb_binary_fields ...

reb_binary_field:
    type: END
    size: 0

reb_simulationarchive_blob:
    index: 1
    offset_prev: 256 (offset to the previous blob)
    offset_next: 256 (offset to the next blob)

reb_binary_field:
    type: DT
    size: 8 bytes

8 bytes of data representing the value of DT

... more reb_binary_fields ...

reb_binary_field:
    type: END
    size: 0

reb_simulationarchive_blob:
    index: 2
    offset_prev: 256 (offset to the previous blob)
    offset_next: 0 
```

The offsets are also used as a sort of checksum to detect if a binary file has been corrupted (for example because a user ran out of disk space). 
If a binary file is corrupted, REBOUND attempts some magic and will recover the last snapshot which does not appear corrupted.
You will see a warning message when that happens and should proceed with caution (make a backup!). 


## Binary Field Descriptor

REBOUND maintains a list of fields it needs to input/output in order to restore a simulation. 
This list is of type `struct reb_binary_field_descriptor[]` and defined in `output.c` as `reb_binary_field_descriptor_list`.
A single struct `reb_binary_field_descriptor` contains the information to input/output one REBOUND field, for example the current simulation time `t`:

```c
    struct reb_binary_field_descriptor fd_t = { 0, REB_DOUBLE, "t", offsetof(struct reb_simulation, t), 0, 0};
```
The first number is a unique identifier (in this case 0). The second entry is the type of data, in this case a single double precision floating point number. The third entry is a string used to identify the field. This is only used when generating human-readable output and is typically the same as the variable name in C. The next entry is the offset of where this variable is stored relative to the beginning of the simulation structure. 

REBOUND also supports array like fields. For example consider the `particles` field:
```c
    struct reb_binary_field_descriptor fd_particles = { 85, REB_POINTER, "particles", offsetof(struct reb_simulation, particles), offsetof(struct reb_simulation, N), sizeof(struct reb_particle)};
```

The second to last entry lists the offset of the a variable in the `reb_simulation` structure that determines the number of array elements. In this case the number of particles. The last entry is the size of a single element. In this case, the size of one `reb_particle`.

If you add an additional field to the `reb_simulation` struct and you want to write it to a binary file and read it back in, then you need to add an entry to `reb_binary_field_descriptor_list`.

