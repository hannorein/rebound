# Binary Format

REBOUND comes with its one binary format.
The binary format allows you to store a current simulation state to a file.
The Simulation Archive is an extension of that format which allows you to store multiple snapshots in one file.
This page explains the details of the binary format.
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
    int16_t offset_prev;
    int16_t offset_next;
};
```

## Binary file (one snapshot)
You create a binary file if you save a simulation
=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... setup simulation ...
    reb_output_binary(r, "snapshot.bin");
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    // ... setup simulation ...
    sim.save("snapshot.bin")
    ```
Such a binary file with one snapshot is simply a set of `reb_binaryfield`s followed by one `reb_simulationarchive_blob` at the end:

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
The type is an integer defined in the `rebound.h` file as `enum REB_BINARY_FIELD_TYPE`.
The last binary field is of type `END` to indicate that the snapshot ends here. 


## SimulationArchive file (multiple snapshots)
The binary file above can also be interpeted as a SimulationArchive with one snapshot. 
You can append many (millions!) of snapshots to a binary file.
REBOUND only stores data that has changed since the original snapshot (typically the particle data, time, etc).
This allows for a very compact file size, while still maintaining bit-wise reproducibility. 

Each snapshot is separated by a `reb_simulationarchive_blob`. 
The blob contains the offset to the previous and next blobs. 
This allows REBOUND to quickly jump from one blob in the archive to the next.
Between the blobs are the same `reb_binary_field`s we already encountered for a binary file with one snapshot.
Thus, a SimulationArchive file with multiple snapshots looks something like this:

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
