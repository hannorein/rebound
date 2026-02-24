/**
 * @file 	simulationarchive.h
 * @brief 	Tools for creating and reading a Simulationarchive binary file.
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
enum REB_BINARYDATA_ERROR_CODE; // Forward declaration.

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

// Used in the binary file to identify data blobs (snapshots)
struct reb_simulationarchive_blob {  
    int32_t index;                   // Index of previous blob (binary file is 0, first blob is 1)
    int32_t offset_prev;             // Offset to beginning of previous blob (size of previous blob).
    int32_t offset_next;             // Offset to end of following blob (size of following blob).
};


// Used by python.
DLLEXPORT void reb_simulation_create_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, int64_t snapshot, enum REB_BINARYDATA_ERROR_CODE* warnings);
DLLEXPORT void reb_simulationarchive_init_from_buffer_with_messages(struct reb_simulationarchive* sa, char* buf, size_t size, struct reb_simulationarchive* sa_index, enum REB_BINARYDATA_ERROR_CODE* warnings);
DLLEXPORT void reb_simulationarchive_create_from_file_with_messages(struct reb_simulationarchive* sa, const char* filename,  struct reb_simulationarchive* sa_index, enum REB_BINARYDATA_ERROR_CODE* warnings);
DLLEXPORT void reb_simulationarchive_free_pointers(struct reb_simulationarchive* sa);

void reb_simulationarchive_heartbeat(struct reb_simulation* const r);  ///< Internal function to handle outputs for the Simulationarchive.
void reb_simulationarchive_create_from_file_with_messages(struct reb_simulationarchive* sa, const char* filename, struct reb_simulationarchive* sa_shape, enum REB_BINARYDATA_ERROR_CODE* warnings); ///< Internal function to read one snapshot from a simulationarchive.


#endif 	// SIMULATIONARCHIVE_H
