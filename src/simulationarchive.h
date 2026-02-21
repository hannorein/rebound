/**
 * @file 	simulationarchive.h
 * @brief 	Tools for creating and readin a Simulationarchive binary file.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2016 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef SIMULATIONARCHIVE_H
#define SIMULATIONARCHIVE_H

#include <stdint.h>

// Simulationarchive structure
struct reb_simulationarchive{
    FILE* inf;                      // File pointer (will be kept open)
    char* filename;                 // Filename of open file. This is NULL if this is a memory-mapped file (using fmemopen)
    int version;                    // Simulationarchive version
    int reb_version_major;          // Major REBOUND Version used to save SA
    int reb_version_minor;          // Minor REBOUND Version used to save SA
    int reb_version_patch;          // Patch REBOUND Version used to save SA
    double auto_interval;           // Interval setting used to create SA (if used)
    double auto_walltime;           // Walltime setting used to create SA (if used)
    uint64_t auto_step;             // Steps in-between SA snapshots (if used)
    int64_t nblobs;                 // Total number of snapshots (including initial binary)
    uint64_t* offset;               // Index of offsets in file (length nblobs)
    double* t;                      // Index of simulation times in file (length nblobs)
};

// Used in the binary file to identify data blobs
struct reb_simulationarchive_blob {  
    int32_t index;                   // Index of previous blob (binary file is 0, first blob is 1)
    int32_t offset_prev;             // Offset to beginning of previous blob (size of previous blob).
    int32_t offset_next;             // Offset to end of following blob (size of following blob).
};

enum REB_BINARY_FIELD_DTYPE {
    REB_DOUBLE = 0,
    REB_INT = 1,
    REB_UINT = 2,                // Same as UINT32
    REB_UINT32 = 3,
    REB_INT64 = 4,
    REB_UINT64 = 5,
    // REB_ULONGLONG = 6,        // No longer used. Using explicit lengths instead.
    REB_VEC3D = 7,
    REB_PARTICLE = 8,
    REB_POINTER = 9,
    REB_POINTER_ALIGNED = 10,    // memory aligned to 64 bit boundary for AVX512
    REB_DP7 = 11,                // Special datatype for IAS15
    REB_OTHER = 12,              // Fields that need special treatment during input and/or output
    REB_FIELD_END = 13,          // Special type to indicate end of blob
    REB_FIELD_NOT_FOUND = 14,    // Special type used to throw error messages
    REB_PARTICLE4 = 15,          // Used for WHFast512
    REB_POINTER_FIXED_SIZE = 16, // A pointer with a fixed size.
    REB_CHARP_LIST = 17,         // A list of NULL terminated strings (char**).
};

// Binary field descriptors are used to identify data blobs in simulationarchives.
struct reb_binarydata_field_descriptor {
    uint32_t type;          // Unique id for each field. Should not change between versions. Ids should not be reused.
    enum REB_BINARY_FIELD_DTYPE dtype; // Datatype (note: not the same as type)
    char name[1024];
    size_t offset;              // Offset of the storage location relative to the beginning of reb_simulation
    size_t offset_N;            // Offset of the storage location for the size relative to the beginning of reb_simulation
    size_t element_size;        // Size in bytes of each element (only used for pointers, dp7, etc).
};

struct reb_binary_field { // This structure is used to save and load binary files.
    uint32_t type;  // type as given by reb_binarydata_field_descriptor
    uint64_t size;  // Size in bytes of field (only counting what follows, not the binary field, itself).
};

// Used by python for testing.
DLLEXPORT extern const struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_list[]; // List of blobs. Implemented in output.c
DLLEXPORT struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_type(int type);
DLLEXPORT struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_name(const char* name);

void reb_simulationarchive_heartbeat(struct reb_simulation* const r);  ///< Internal function to handle outputs for the Simulationarchive.
void reb_simulationarchive_create_from_file_with_messages(struct reb_simulationarchive* sa, const char* filename, struct reb_simulationarchive* sa_shape, enum reb_simulation_binary_error_codes* warnings); ///< Internal function to read one snapshot from a simulationarchive.


#endif 	// SIMULATIONARCHIVE_H
