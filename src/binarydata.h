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
    REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR = 2048,
};

enum REB_BINARYDATA_OUTPUT {
    REB_BINARYDATA_OUTPUT_NONE = 0,
    REB_BINARYDATA_OUTPUT_PRINT = 1,
    REB_BINARYDATA_OUTPUT_STREAM = 2,
    REB_BINARYDATA_OUTPUT_BUFFER = 3
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
int reb_binarydata_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, enum REB_BINARYDATA_OUTPUT output_option);

// Read all fields from inf stream into r. 
void reb_binarydata_input_fields(struct reb_simulation* r, FILE* inf, enum REB_BINARYDATA_ERROR_CODE* warnings); 

// Process any errors that might have occurred while reading binary data.
struct reb_simulation* reb_binarydata_process_warnings(struct reb_simulation* r, enum REB_BINARYDATA_ERROR_CODE warnings);

// Write the simulation to a memory buffer (simulationarchive format).
REB_API void reb_binarydata_simulation_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep);

#define REB_AS_STRUCT_MEMBER(prefix, value, name) {value, #name}, 
#define REB_GENERATE_ENUM_DESCRIPTORS(LIST) \
    (struct reb_binarydata_enum_descriptor[]){ \
        LIST(REB_AS_STRUCT_MEMBER, LIST)\
        {0} /* Null terminated*/ \
    }

// This structure is written/read to files. Precedes the actual data. 
struct reb_binarydata_field { 
    uint32_t type;  // type as given by reb_binarydata_field_descriptor
    uint64_t size;  // Size in bytes of field (only counting what follows, not including reb_binarydata_field itself).
};

// List of all possible input/ouput fields
REB_API extern const struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_list[];
// Helper functions to find descriptor data. Used by python.
REB_API struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_type(uint32_t type);
REB_API struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_name(const char* name);

#endif // _BINARYDATA_H
