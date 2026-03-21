/**
 * @file 	binarydata.h
 * @brief 	Routines for output, input and comparison of simulations in binary format.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2026 Hanno Rein
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
#ifndef _BINARYDATA_H
#define _BINARYDATA_H
#include <stdio.h>

// Possible datatypes for reb_binarydata_field.
enum REB_BINARYDATA_DTYPE {
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
    REB_SIZE_T = 18,
};

// Possible errors that might occur during binary file reading.
enum REB_BINARYDATA_ERROR_CODE {
    REB_BINARYDATA_WARNING_NONE = 0,
    REB_BINARYDATA_ERROR_NOFILE = 1,
    REB_BINARYDATA_WARNING_VERSION = 2,
    REB_BINARYDATA_WARNING_POINTERS = 4,
    REB_BINARYDATA_WARNING_PARTICLES = 8,
    REB_BINARYDATA_ERROR_FILENOTOPEN = 16,
    REB_BINARYDATA_ERROR_OUTOFRANGE = 32,
    REB_BINARYDATA_ERROR_SEEK = 64,
    REB_BINARYDATA_WARNING_FIELD_UNKNOWN = 128,
    REB_BINARYDATA_ERROR_INTEGRATOR = 256,
    REB_BINARYDATA_WARNING_CORRUPTFILE = 512,
    REB_BINARYDATA_ERROR_OLD = 1024,
};

// Compares two simulations, stores difference in buffer.
//
// output_option:
// - If set to 0, differences are written to bufp in the form of reb_binarydata_field structs. 
// - If set to 1, differences are printed on the screen. 
// - If set to 2, only the return value indicates any differences.
// - If set to 3, differences are written to bufp in a human readable form.
//
// returns value:  0 is returned if the simulations do not differ (are equal). 1 is return if they differ.
int reb_binarydata_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option);

// Read all fields from inf stream into r. 
void reb_binarydata_input_fields(struct reb_simulation* r, FILE* inf, enum REB_BINARYDATA_ERROR_CODE* warnings); 

// Process any errors that might have occurred while reading binary data.
struct reb_simulation* reb_binarydata_process_warnings(struct reb_simulation* r, enum REB_BINARYDATA_ERROR_CODE warnings);

// Write the simulation to a memory buffer (simulationarchive format).
DLLEXPORT void reb_binarydata_simulation_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep);

// Binary field descriptors are used to identify data blobs in simulationarchives.
struct reb_binarydata_field_descriptor {
    uint32_t type;              // Unique id for each field. Should not change between versions. Ids should not be reused.
    enum REB_BINARYDATA_DTYPE dtype; // Datatype (note: not the same as type)
    char name[1024];            // Null terminated human readable name.
    size_t offset;              // Offset of the storage location relative to the beginning of reb_simulation
    size_t offset_N;            // Offset of the storage location for the size relative to the beginning of reb_simulation
    size_t element_size;        // Size in bytes of each element (only used for pointers, dp7, etc).
};

// This structure is written/read to files. Precedes the actual data. 
struct reb_binarydata_field { 
    uint32_t type;  // type as given by reb_binarydata_field_descriptor
    uint64_t size;  // Size in bytes of field (only counting what follows, not including reb_binarydata_field itself).
};

// List of all possible input/ouput fields
DLLEXPORT extern const struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_list[];
// Helper functions to find descriptor data. Used by python.
DLLEXPORT struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_type(uint32_t type);
DLLEXPORT struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_name(const char* name);

#endif // _BINARYDATA_H
