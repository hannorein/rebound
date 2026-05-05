/**
 * @file    binarydata.c
 * @brief 	Routines for output, input and comparison of simulations in binary format.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
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
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <inttypes.h>
#include "rebound.h"
#include "rebound_internal.h"
#include "particle.h"
#include "tools.h"
#include "tree.h"
#include "output.h"
#include "binarydata.h"
#include "simulationarchive.h"

const uint64_t reb_binarydata_header = 0x20444E554F424552; // Corresponds to the first few ASCII characters in binary file

// Null terminated list of REBOUND parameters to be written to a file.
// Modify this list if you wish to input/output additional fields in the reb_simulation structure.
const struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_list[]= {
    { REB_DOUBLE,       "t",                            offsetof(struct reb_simulation, t), 0, 0, 0}, // used to be id 0
    { REB_DOUBLE,       "G",                            offsetof(struct reb_simulation, G), 0, 0, 0},
    { REB_DOUBLE,       "softening",                    offsetof(struct reb_simulation, softening), 0, 0, 0},
    { REB_DOUBLE,       "dt",                           offsetof(struct reb_simulation, dt), 0, 0, 0},
    { REB_SIZE_T,       "N",                            offsetof(struct reb_simulation, N), 0, 0, 0},
    { REB_SIZE_T,       "N_var",                        offsetof(struct reb_simulation, N_var), 0, 0, 0},
    { REB_SIZE_T,       "N_active",                     offsetof(struct reb_simulation, N_active), 0, 0, 0},
    { REB_INT,          "testparticle_type",            offsetof(struct reb_simulation, testparticle_type), 0, 0, 0},
    { REB_DOUBLE,       "opening_angle2",               offsetof(struct reb_simulation, opening_angle2), 0, 0, 0},
    { REB_INT,          "status",                       offsetof(struct reb_simulation, status), 0, 0, 0},
    { REB_INT,          "exact_finish_time",            offsetof(struct reb_simulation, exact_finish_time), 0, 0, 0},
    { REB_UINT,         "force_is_velocity_dependent",  offsetof(struct reb_simulation, force_is_velocity_dependent), 0, 0, 0},
    { REB_UINT,         "gravity_ignore_terms",         offsetof(struct reb_simulation, gravity_ignore_terms), 0, 0, 0},
    { REB_DOUBLE,       "output_timing_last",           offsetof(struct reb_simulation, output_timing_last), 0, 0, 0},
    { REB_INT,          "save_messages",                offsetof(struct reb_simulation, save_messages), 0, 0, 0},
    { REB_DOUBLE,       "exit_max_distance",            offsetof(struct reb_simulation, exit_max_distance), 0, 0, 0},
    { REB_DOUBLE,       "exit_min_distance",            offsetof(struct reb_simulation, exit_min_distance), 0, 0, 0},
    { REB_DOUBLE,       "usleep",                       offsetof(struct reb_simulation, usleep), 0, 0, 0},
    { REB_INT,          "track_energy_offset",          offsetof(struct reb_simulation, track_energy_offset), 0, 0, 0},
    { REB_DOUBLE,       "energy_offset",                offsetof(struct reb_simulation, energy_offset), 0, 0, 0},
    { REB_DOUBLE,       "root_size",                    offsetof(struct reb_simulation, root_size), 0, 0, 0},
    { REB_SIZE_T,       "N_root_x",                     offsetof(struct reb_simulation, N_root_x), 0, 0, 0},
    { REB_SIZE_T,       "N_root_y",                     offsetof(struct reb_simulation, N_root_y), 0, 0, 0},
    { REB_SIZE_T,       "N_root_z",                     offsetof(struct reb_simulation, N_root_z), 0, 0, 0},
    { REB_INT,          "N_ghost_x",                    offsetof(struct reb_simulation, N_ghost_x), 0, 0, 0},
    { REB_INT,          "N_ghost_y",                    offsetof(struct reb_simulation, N_ghost_y), 0, 0, 0},
    { REB_INT,          "N_ghost_z",                    offsetof(struct reb_simulation, N_ghost_z), 0, 0, 0},
    { REB_DOUBLE,       "minimum_collision_velocity",   offsetof(struct reb_simulation, minimum_collision_velocity), 0, 0, 0},
    { REB_DOUBLE,       "collisions_plog",              offsetof(struct reb_simulation, collisions_plog), 0, 0, 0},
    { REB_INT64,        "collisions_log_n",             offsetof(struct reb_simulation, collisions_log_n), 0, 0, 0},
    { REB_INT,          "calculate_megno",              offsetof(struct reb_simulation, calculate_megno), 0, 0, 0},
    { REB_DOUBLE,       "megno_Ys",                     offsetof(struct reb_simulation, megno_Ys), 0, 0, 0},
    { REB_DOUBLE,       "megno_Yss",                    offsetof(struct reb_simulation, megno_Yss), 0, 0, 0},
    { REB_DOUBLE,       "megno_cov_Yt",                 offsetof(struct reb_simulation, megno_cov_Yt), 0, 0, 0},
    { REB_DOUBLE,       "megno_var_t",                  offsetof(struct reb_simulation, megno_var_t), 0, 0, 0},
    { REB_DOUBLE,       "megno_mean_t",                 offsetof(struct reb_simulation, megno_mean_t), 0, 0, 0},
    { REB_DOUBLE,       "megno_mean_Y",                 offsetof(struct reb_simulation, megno_mean_Y), 0, 0, 0},
    { REB_DOUBLE,       "megno_initial_t",              offsetof(struct reb_simulation, megno_initial_t), 0, 0, 0},
    { REB_INT64,        "megno_n",                      offsetof(struct reb_simulation, megno_n), 0, 0, 0},
    { REB_DOUBLE,       "simulationarchive_auto_interval", offsetof(struct reb_simulation, simulationarchive_auto_interval), 0, 0, 0},
    { REB_DOUBLE,       "simulationarchive_auto_walltime", offsetof(struct reb_simulation, simulationarchive_auto_walltime), 0, 0, 0},
    { REB_DOUBLE,       "simulationarchive_next",       offsetof(struct reb_simulation, simulationarchive_next), 0, 0, 0},
    { REB_INT,          "collision",                    offsetof(struct reb_simulation, collision), 0, 0, 0},
    { REB_STRING,       "integrator.name",              offsetof(struct reb_simulation, integrator.name), 0, 0, 0},
    { REB_INT,          "boundary",                     offsetof(struct reb_simulation, boundary), 0, 0, 0},
    { REB_INT,          "gravity",                      offsetof(struct reb_simulation, gravity), 0, 0, 0},
    { REB_DOUBLE,       "OMEGA",                        offsetof(struct reb_simulation, OMEGA), 0, 0, 0},
    { REB_DOUBLE,       "OMEGAZ",                       offsetof(struct reb_simulation, OMEGAZ), 0, 0, 0},
    { REB_UINT,         "is_synchronized",              offsetof(struct reb_simulation, is_synchronized), 0, 0, 0},
    { REB_UINT,         "did_modify_particles",         offsetof(struct reb_simulation, did_modify_particles), 0, 0, 0},
    { REB_POINTER,      "particles",                    offsetof(struct reb_simulation, particles), offsetof(struct reb_simulation, N), sizeof(struct reb_particle), 0},
    { REB_POINTER,      "particles_var",                offsetof(struct reb_simulation, particles_var), offsetof(struct reb_simulation, N_var), sizeof(struct reb_particle), 0},
    { REB_POINTER,      "var_config",                   offsetof(struct reb_simulation, var_config), offsetof(struct reb_simulation, N_var_config), sizeof(struct reb_variational_configuration), 0},
    { REB_INT,          "simulationarchive_version",    offsetof(struct reb_simulation, simulationarchive_version), 0, 0, 0},
    { REB_DOUBLE,       "walltime",                     offsetof(struct reb_simulation, walltime), 0, 0, 0},
    { REB_DOUBLE,       "walltime_last_steps",          offsetof(struct reb_simulation, walltime_last_steps), 0, 0, 0},
    { REB_UINT32,       "python_unit_l",                offsetof(struct reb_simulation, python_unit_l), 0, 0, 0},
    { REB_UINT32,       "python_unit_m",                offsetof(struct reb_simulation, python_unit_m), 0, 0, 0},
    { REB_UINT32,       "python_unit_t",                offsetof(struct reb_simulation, python_unit_t), 0, 0, 0},
    { REB_UINT64,       "simulationarchive_auto_step",  offsetof(struct reb_simulation, simulationarchive_auto_step), 0, 0, 0},
    { REB_UINT64,       "simulationarchive_next_step",  offsetof(struct reb_simulation, simulationarchive_next_step), 0, 0, 0},
    { REB_UINT64,       "steps_done",                   offsetof(struct reb_simulation, steps_done), 0, 0, 0},
    { REB_DOUBLE,       "dt_last_done",                 offsetof(struct reb_simulation, dt_last_done), 0, 0, 0},
    { REB_UINT,         "rand_seed",                    offsetof(struct reb_simulation, rand_seed), 0, 0, 0},
    { REB_INT,          "testparticle_hidewarnings",    offsetof(struct reb_simulation, testparticle_hidewarnings), 0, 0, 0},
    { REB_POINTER,      "display_settings",             offsetof(struct reb_simulation, display_settings), SIZE_MAX, sizeof(struct reb_display_settings), 0},  // Note: SIZE_MAX means 1 element if pointer not NULL
    { REB_CHARP_LIST,   "name_list",                    offsetof(struct reb_simulation, name_list), offsetof(struct reb_simulation, N_name_list), 0, 0},
    { REB_OTHER,        "functionpointers", 0, 0, 0, 0},
    { REB_OTHER,        "sablob", 0, 0, 0, 0},
    { REB_OTHER,        "end", 0, 0, 0, 0},
    { REB_OTHER,        "header", 0, 0, 0, 0},
    {0} // Null terminated.
};


// This is a custom implementation of a dynamic memory buffer stream. 
// This is used as a replacement for open_memstream which is not portable.
static void write_to_stream(char** bufp, size_t* allocatedsize, size_t* sizep, void* restrict data, size_t size){
    // Increase size
    int increased = 0;
    while (*allocatedsize==0 || (*sizep)+size>(*allocatedsize)){
        increased = 1;
        *allocatedsize = (*allocatedsize) ? (*allocatedsize) * 2 : 32;
    }
    if (increased){
        *bufp = realloc(*bufp,*allocatedsize);
    }
    // Copy data to buffer
    memcpy((*bufp)+(*sizep),data,size);
    *sizep += size;
}

static const struct reb_binarydata_field_descriptor* reb_binarydata_field_descriptor_for_name_in_list(const struct reb_binarydata_field_descriptor* list, const char* name){
    if (!list) return NULL;
    for(size_t i=0; list[i].name[0]; i++){
        if (strcmp(list[i].name, name)==0){
            return &list[i];
        }
    }
    return NULL;
}

// Returns a field descriptor with matching name
// Modifies field descriptor so that name is full qualifying name including prefix.
// Modifies field descriptor so that offset is actual memory address if r is given. 
struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_name(const struct reb_simulation * const r, const char* name){
    const struct reb_binarydata_field_descriptor* fd = NULL;
    char* name_sub;
    // Check if this is an integrator field.
    if (strncmp("integrator.", name, 11)==0 && (name_sub = strchr(name+11,'.'))){
        // First, check current integrator.
        if (r && r->integrator.callbacks.field_descriptor_list){
            fd = reb_binarydata_field_descriptor_for_name_in_list(r->integrator.callbacks.field_descriptor_list, name_sub+1);
        }else{
            // Look through all built-in integrators
#define X(iname) if (!fd) {fd = reb_binarydata_field_descriptor_for_name_in_list(reb_integrator_##iname.field_descriptor_list, name_sub+1);}
            REB_AVAILABLE_INTEGRATORS
#undef X
                // Look through all custom integrators
                size_t Nc = 0;
                if (reb_integrator_configurations_custom){
                    while(reb_integrator_configurations_custom[Nc].name){
                        fd = reb_binarydata_field_descriptor_for_name_in_list(reb_integrator_configurations_custom[Nc].callbacks.field_descriptor_list, name_sub+1);
                        if (fd) break; // found field
                        Nc++;
                    }
                }
        }
        if (fd){
            struct reb_binarydata_field_descriptor fd_integrator = *fd;
            if (r && r->integrator.state){
                fd_integrator.offset += (size_t)r->integrator.state;
                fd_integrator.offset_N += (size_t)r->integrator.state;
            }
            strcpy(fd_integrator.name, name);
            return fd_integrator;
        } 
    }
    
    fd = reb_binarydata_field_descriptor_for_name_in_list(reb_binarydata_field_descriptor_list, name); 
    if (fd){
        struct reb_binarydata_field_descriptor fd_simulation = *fd;
        if (r){
            fd_simulation.offset += (size_t)r;
            fd_simulation.offset_N += (size_t)r;
        }
        return fd_simulation;
    }
            printf("here 2\n");
    reb_simulation_error((struct reb_simulation*)r, "Could not find field descriptor for name.");
    struct reb_binarydata_field_descriptor bfd = {
        .dtype = REB_FIELD_NOT_FOUND,
        .name = "field_not_found",
    };
    return bfd;
}

// Like asprintf, but append to bufp
static void asprintf_append_to_bufp(char** bufp, const char* format, ...){
    char* buf = NULL;
    va_list args;
    va_start(args, format);
    vasprintf(&buf, format, args);
    va_end(args);
    if (bufp){
        *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
        strcat(*bufp,buf);
    }else{
        printf("%s",buf);
    }
    free(buf);
}

// Helper function to print out binary data in human readable form.
static void asprintf_append_to_bufp_type(char** bufp, enum REB_BINARYDATA_DTYPE dtype, char* pointer, size_t dsize){
    switch (dtype){
        case REB_DOUBLE:
            asprintf_append_to_bufp(bufp, "%e",*(double*)(pointer));
            break;
        case REB_INT:
            asprintf_append_to_bufp(bufp, "%d",*(int*)(pointer));
            break;
        case REB_SIZE_T:
            asprintf_append_to_bufp(bufp, "%zu",*(size_t*)(pointer));
            break;
        case REB_UINT:
            asprintf_append_to_bufp(bufp, "%u",*(unsigned int*)(pointer));
            break;
        case REB_UINT32:
            asprintf_append_to_bufp(bufp, "%" PRIu32,*(uint32_t*)(pointer)); // PRIu32 defined in inttypes.h
            break;
        case REB_INT64:
            asprintf_append_to_bufp(bufp, "%" PRId64,*(int64_t*)(pointer));
            break;
        case REB_UINT64:
            asprintf_append_to_bufp(bufp, "%" PRIu64,*(uint64_t*)(pointer));
            break;
        case REB_STRING:
            if (pointer && *pointer){
                asprintf_append_to_bufp(bufp, "\"%s\"",pointer);
            }else{
                asprintf_append_to_bufp(bufp, "NULL");
            }
            break;
        default:
            asprintf_append_to_bufp(bufp, "(%zu bytes, values not printed)", dsize);
            break;
    }
}

// Compares two simulations in buffers.
// Returns 0 if the buffers contain the same simulation data. 
// Supports different output options.
int reb_binarydata_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, enum REB_BINARYDATA_OUTPUT output_option){
    if (!buf1 || !buf2 || size1<64 || size2<64){
        printf("Cannot read input buffers.\n");
        return 0;
    }

    int are_different = 0;

    if (output_option==REB_BINARYDATA_OUTPUT_STREAM){
        *bufp = NULL;
        *sizep = 0;
    }
    if (output_option==REB_BINARYDATA_OUTPUT_BUFFER){
        *bufp = malloc(sizeof(char));
        *bufp[0] = '\0';
    }
    size_t allocatedsize = 0;

    // Header.
    if(memcmp(buf1,buf2,64)!=0 && output_option==REB_BINARYDATA_OUTPUT_PRINT){
        printf("Header in binary files are different.\n");
    }

    size_t pos1 = 64;
    size_t pos2 = 64;
    struct reb_binarydata_field field1;
    struct reb_binarydata_field field2;
    char* name1;
    char* name2;

    while(1){
        if (pos1+sizeof(struct reb_binarydata_field)>size1) break;
        memcpy(&field1, buf1+pos1, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos1 += sizeof(struct reb_binarydata_field);
        name1 = buf1+pos1;
        pos1 += field1.size_name;
        if (strcmp(name1, "end")==0){
            break;
        }
        if (pos2+sizeof(struct reb_binarydata_field)>size2) pos2 = 64;
        memcpy(&field2, buf2+pos2, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos2 += sizeof(struct reb_binarydata_field);
        name2 = buf2+pos2;
        pos2 += field2.size_name;

        // Fields might not be in the same order.
        if (strcmp(name1, name2)){
            // Will search for element in buf2, starting at beginning just past header
            // Note that we ignore all ADDITIONAL fields in buf2 that were not present in buf1 
            pos2 = 64;
            int notfound = 0; 
            while(1) {
                if (pos2+sizeof(struct reb_binarydata_field)>size2){
                    notfound = 1;
                    break;
                }
                memcpy(&field2, buf2+pos2, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
                pos2 += sizeof(struct reb_binarydata_field);
                name2 = buf2+pos2;
                pos2 += field2.size_name;
                if(strcmp(name2, "end")==0){
                    notfound = 1;
                    break;
                }
                if (strcmp(name1, name2)==0){
                    break; // found!!
                }else{
                    pos2 += field2.size_data; //skip
                }
            };
            if (notfound == 1){
                are_different = 1.;
                switch (output_option){
                    case REB_BINARYDATA_OUTPUT_STREAM:
                        write_to_stream(bufp, &allocatedsize, sizep, &field1,sizeof(struct reb_binarydata_field));
                        write_to_stream(bufp, &allocatedsize, sizep, name1,field1.size_name);
                        break;
                    case REB_BINARYDATA_OUTPUT_PRINT:
                    case REB_BINARYDATA_OUTPUT_BUFFER:
                        {
                            const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_name(NULL, name1);
                            asprintf_append_to_bufp(bufp, "%s:\n" REB_STR_RED "< ",name1);
                            asprintf_append_to_bufp_type(bufp, fd.dtype, buf1+pos1, field1.size_data);
                            asprintf_append_to_bufp(bufp, REB_STR_RESET "\n");
                        }
                        break;
                    case REB_BINARYDATA_OUTPUT_NONE:
                        break;
                }           
                // Set offsets for next search
                pos2 = 64;
                pos1 += field1.size_data;
                continue;
            }
        }
        // Can assume field1 and field2 have the same name from here on
        if (pos1+field1.size_data>size1) printf("Corrupt binary file buf1.\n");
        if (pos2+field2.size_data>size2) printf("Corrupt binary file buf2.\n");
        int fields_differ = 0;
        if (field1.size_data==field2.size_data){
            if (strcmp(name1, "particles")==0){ // TODO Should be dtype== REB_PARTICLE so it works for all particle arrays
                struct reb_particle* pb1 = (struct reb_particle*)(buf1+pos1);
                struct reb_particle* pb2 = (struct reb_particle*)(buf2+pos2);
                for (size_t i=0;i<field1.size_data/sizeof(struct reb_particle);i++){
                    struct reb_particle p1;
                    struct reb_particle p2;
                    memcpy(&p1, pb1+i, sizeof(struct reb_particle)); // need copy because of 8 byte alignment requirement
                    memcpy(&p2, pb2+i, sizeof(struct reb_particle)); // need copy because of 8 byte alignment requirement
                    fields_differ |= reb_particle_cmp(p1,p2);
                }
            }else{
                if (memcmp(buf1+pos1,buf2+pos2,field1.size_data)!=0){
                    fields_differ = 1;
                }
            }
        }else{
            fields_differ = 1;
        }
        if(fields_differ){
            if (strcmp(name1, "walltime")!=0){
                // Ignore the walltime fields, but only for the return value (print it out)
                // Typically we do not care about this field when comparing simulations.
                are_different = 1.;
            }
            switch (output_option){
                case REB_BINARYDATA_OUTPUT_STREAM:
                    write_to_stream(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binarydata_field));
                    write_to_stream(bufp, &allocatedsize, sizep, name2,field2.size_name);
                    write_to_stream(bufp, &allocatedsize, sizep, buf2+pos2,field2.size_data);
                    break;
                case REB_BINARYDATA_OUTPUT_PRINT:
                case REB_BINARYDATA_OUTPUT_BUFFER:
                    {
                        const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_name(NULL, name1);
                        asprintf_append_to_bufp(bufp, "%s:\n" REB_STR_RED "< ",name1);
                        asprintf_append_to_bufp_type(bufp, fd.dtype, buf1+pos1, field1.size_data);
                        asprintf_append_to_bufp(bufp, REB_STR_RESET "\n---\n" REB_STR_GREEN "> ");
                        asprintf_append_to_bufp_type(bufp, fd.dtype, buf2+pos2, field2.size_data);
                        asprintf_append_to_bufp(bufp, REB_STR_RESET "\n");
                    }
                    break;
                case REB_BINARYDATA_OUTPUT_NONE:
                    break;
            }
        }
        pos1 += field1.size_data;
        pos2 += field2.size_data;
    }
    // Search for fields which are present in buf2 but not in buf1
    pos1 = 64;
    pos2 = 64;
    while(1){
        if (pos2+sizeof(struct reb_binarydata_field)>size2) break;
        memcpy(&field2, buf2+pos2, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos2 += sizeof(struct reb_binarydata_field);
        name2 = buf2+pos2;
        pos2 += field2.size_name;
        if (strcmp(name2, "end")==0){
            break;
        }
        if (pos1+sizeof(struct reb_binarydata_field)>size1) pos1 = 64;
        memcpy(&field1, buf1+pos1, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos1 += sizeof(struct reb_binarydata_field);
        name1 = buf1+pos1;
        pos1 += field1.size_name;

        if (strcmp(name1, name2)==0){
            // Not a new field. Skip.
            pos1 += field1.size_data;
            pos2 += field2.size_data;
            continue;
        }
        // Fields might not be in the same order.
        // Will search for element in buf1, starting at beginning just past header
        pos1 = 64;
        int notfound = 0; 
        while(1) {
            if (pos1+sizeof(struct reb_binarydata_field)>size1){
                notfound = 1;
                break;
            }
            memcpy(&field1, buf1+pos1, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
            pos1 += sizeof(struct reb_binarydata_field);
            name1 = buf1+pos1;
            pos1 += field1.size_name;
            if(strcmp(name1, "end")==0){
                notfound = 1;
                break;
            }
            if (strcmp(name1, name2)==0){
                break; // found it, not new
            }else{
                // not found, try next
                pos1 += field1.size_data;
            }
        };
        if (notfound == 0){
            // Not a new field. Skip.
            pos1 = 64;
            pos2 += field2.size_data;
            continue;
        }

        are_different = 1.;
        switch (output_option){
            case REB_BINARYDATA_OUTPUT_STREAM:
                write_to_stream(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, &allocatedsize, sizep, name2,field2.size_name);
                write_to_stream(bufp, &allocatedsize, sizep, buf2+pos2,field2.size_data);
                break;
            case REB_BINARYDATA_OUTPUT_PRINT:
            case REB_BINARYDATA_OUTPUT_BUFFER:
                {
                    const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_name(NULL, name2);
                    asprintf_append_to_bufp(bufp, "%s:\n" REB_STR_GREEN "> ",name2);
                    asprintf_append_to_bufp_type(bufp, fd.dtype, buf2+pos2, field2.size_data);
                    asprintf_append_to_bufp(bufp, REB_STR_RESET "\n");
                }
                break;
            case REB_BINARYDATA_OUTPUT_NONE:
                break;
        }
        pos1 = 64;
        pos2 += field2.size_data;
    }

    return are_different;
}

// Output all fields from one field_descriptor list
static void output_fields_from_list(char** bufp, size_t* current_pos, size_t* allocatedsize, const struct reb_binarydata_field_descriptor* fd_list, char* base_address, const char* prefix){
    if (!fd_list) return; 
    char name[1024];

    for (size_t i=0; fd_list[i].name[0]; i++){
        struct reb_binarydata_field_descriptor fd = fd_list[i];
        if (prefix && prefix[0]){
            strcpy(name, prefix);
            strcat(name, ".");
            strcat(name, fd.name);
        }else{
            strcpy(name, fd.name);
        }
        size_t size_name = strlen(name)+1;
        struct reb_binarydata_field field = {.size_name = size_name};

        char** pointer_to_value = (char**)(base_address + fd.offset);
        switch (fd.dtype){
            case REB_DOUBLE: 
                field.size_data = sizeof(double);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_INT: 
                field.size_data = sizeof(int);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_SIZE_T: 
                field.size_data = sizeof(size_t);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_UINT: 
                field.size_data = sizeof(unsigned int);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_UINT32: 
                field.size_data = sizeof(uint32_t);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_INT64:
                field.size_data = sizeof(int64_t);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_UINT64:
                field.size_data = sizeof(uint64_t);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_VEC3D:
                field.size_data = sizeof(struct reb_vec3d);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_PARTICLE:
                field.size_data = sizeof(struct reb_particle);
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_STRING: 
                if (*pointer_to_value){
                    field.size_data = strlen(*pointer_to_value)+1;
                    pointer_to_value = (char**)(*pointer_to_value);
                }else{
                    field.size_data = 0;
                }
                write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                write_to_stream(bufp, allocatedsize, current_pos, pointer_to_value, field.size_data);
                break;
            case REB_POINTER:
            case REB_POINTER_ALIGNED:
                {
                    size_t pointer_N = 0;
                    if (fd.offset_N!=SIZE_MAX){ // Dynamic N if offset_N is given
                        pointer_N = *(size_t*)(base_address + fd.offset_N);
                    }else{ // Fixed size_data pointer_to_value.
                        if (*pointer_to_value){
                            pointer_N = 1; // Pointer is not NULL, thus store one element.
                        }
                    }
                    field.size_data = pointer_N * fd.element_size;

                    if (field.size_data){
                        write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                        write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                        write_to_stream(bufp, allocatedsize, current_pos, *pointer_to_value, field.size_data);
                    }
                }
                break;
            case REB_CHARP_LIST:
                {
                    size_t N_list = *((size_t*)(base_address + fd.offset_N));
                    char*** list_p = (char***)pointer_to_value; // pointer to a list of strings
                    size_t serialized_size = 0;
                    for (size_t i=0; i<N_list; i++){
                        // character count + NULL character + original pointer
                        serialized_size += strlen((*list_p)[i])+1+sizeof(char*);
                    }
                    field.size_data = sizeof(char)*serialized_size;

                    if (field.size_data){
                        // This pointer arithmetic will fail on 32 bit architectures.
                        write_to_stream(bufp, allocatedsize, current_pos, &field, sizeof(struct reb_binarydata_field));
                        write_to_stream(bufp, allocatedsize, current_pos, name, field.size_name);
                        for (size_t i=0; i<N_list; i++){
                            write_to_stream(bufp, allocatedsize, current_pos, (*list_p)[i], strlen((*list_p)[i])+1);
                            write_to_stream(bufp, allocatedsize, current_pos, &((*list_p)[i]), sizeof(char*));
                        }
                    }
                }
                break;
            case REB_OTHER:
            case REB_FIELD_NOT_FOUND:
                // Special fields are not written to output.
                break;
        }
    }
}

// Serializes a simulation to a buffer
void reb_binarydata_simulation_to_stream(struct reb_simulation* r, char** bufp, size_t* current_pos){
    if (r->simulationarchive_version<5){
        reb_simulation_error(r, "Simulationarchives with version < 5 are no longer supported.\n");
    }
    size_t allocatedsize = 0;
    *bufp = NULL;
    *current_pos = 0;

    // Output header.
    char header[64] = "\0";
    int cwritten = sprintf(header,"REBOUND Binary File. Version: %s",reb_version_str);
    snprintf(header+cwritten+1,64-cwritten-1,"%s",reb_githash_str);
    write_to_stream(bufp, &allocatedsize, current_pos, header,sizeof(char)*64);

    /// Output all fields
    // Main simulation
    output_fields_from_list(bufp, current_pos, &allocatedsize, reb_binarydata_field_descriptor_list, (char*)r, NULL);
    // Integrator
    if (r->integrator.name[0]){
        char prefix[256] = "integrator.";
        strlcat(prefix, r->integrator.name, sizeof(prefix));
        output_fields_from_list(bufp, current_pos, &allocatedsize, r->integrator.callbacks.field_descriptor_list, (char*)r->integrator.state, prefix);
    }

    // Write function pointer warning flag
    int functionpointersused = 0;
    if (r->coefficient_of_restitution ||
            r->collision_resolve ||
            r->additional_forces ||
            r->heartbeat ||
            r->post_timestep_modifications ||
            r->free_particle_ap){
        functionpointersused = 1;
    }

    struct reb_binarydata_field_descriptor fd_fp = reb_binarydata_field_descriptor_for_name(NULL, "functionpointers");
    struct reb_binarydata_field field_functionp;
    memset(&field_functionp,0,sizeof(struct reb_binarydata_field));
    field_functionp.size_name = strlen(fd_fp.name)+1;
    field_functionp.size_data = sizeof(int);
    write_to_stream(bufp, &allocatedsize, current_pos, &field_functionp, sizeof(struct reb_binarydata_field));
    write_to_stream(bufp, &allocatedsize, current_pos, fd_fp.name, field_functionp.size_name);
    write_to_stream(bufp, &allocatedsize, current_pos, &functionpointersused, field_functionp.size_data);

    // Write last field
    struct reb_binarydata_field_descriptor fd_end = reb_binarydata_field_descriptor_for_name(NULL, "end");
    struct reb_binarydata_field end_field = {.size_name = strlen(fd_end.name)+1, .size_data = 0};
    write_to_stream(bufp, &allocatedsize, current_pos, &end_field, sizeof(struct reb_binarydata_field));
    write_to_stream(bufp, &allocatedsize, current_pos, fd_end.name, end_field.size_name);

    struct reb_simulationarchive_blob blob = {0};
    write_to_stream(bufp, &allocatedsize, current_pos, &blob, sizeof(struct reb_simulationarchive_blob));
}

// Read field data into simulation from file or memory buffer.
void reb_binarydata_input_fields(struct reb_simulation* r, FILE* inf, enum REB_BINARYDATA_ERROR_CODE* warnings){
    struct reb_binarydata_field field;
    char name[1024];
next_field:
    // Loop over all fields
    while(1){

        int numread = (int)fread(&field,sizeof(struct reb_binarydata_field),1,inf);
        if (numread<1){
            goto finish_fields; // End of file
        }
        // Is this a real field or the header?
        if (field.size_name == reb_binarydata_header) {
            int64_t objects = 0;
            const size_t bufsize = 64 - sizeof(struct reb_binarydata_field);
            char readbuf[64], curvbuf[64];
            const char* header = "REBOUND Binary File. Version: ";
            sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binarydata_field), reb_version_str);

            objects += fread(readbuf,sizeof(char),bufsize,inf);
            if (objects < 1){
                *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            }else{
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_BINARYDATA_WARNING_VERSION;
                }
            }
            goto next_field;
        }
        // Try to get name of field
        numread = (int)fread(name,field.size_name,1,inf);
        if (numread<1){
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            goto finish_fields; // End of file
        }
        // Fields that require special handling
        if (strcmp(name, "end")==0){
            goto finish_fields; // End of snapshot
        }
        if (strcmp(name, "functionpointers")==0){
            // Warning for when function pointers were used. 
            // No effect on simulation.
            int fpwarn;
            fread(&fpwarn, field.size_data,1,inf);
            if (fpwarn && warnings){
                *warnings |= REB_BINARYDATA_WARNING_POINTERS;
            }
            goto next_field;
        }


        struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_name(r, name) ;
        char** pointer = (char**)fd.offset;
        switch(fd.dtype){
            case REB_DOUBLE:
            case REB_INT:
            case REB_UINT:
            case REB_UINT32:
            case REB_INT64:
            case REB_UINT64:
            case REB_PARTICLE:
            case REB_VEC3D:
            case REB_SIZE_T:
                // Simple datatypes:
                fread(pointer, field.size_data, 1, inf);
                goto next_field;
                break;
            case REB_POINTER:
            case REB_POINTER_ALIGNED:
                // Read a pointer data type. 
                // 1) reallocate memory
                // 2) read data into memory
                // 3) set N_allocated variable
                {
                    if (field.size_data % fd.element_size){
                        reb_simulation_warning(r, "Inconsistent size_data encountered in binary field.");
                    }
                    size_t* pointer_N = (size_t*)fd.offset_N;
                    if (fd.dtype == REB_POINTER_ALIGNED){
                        if (*pointer) free(*pointer);
#if defined(_WIN32) || !defined(AVX512)
                        // WHFast512 not supported on Windows!
                        *pointer = malloc(field.size_data);
#else 
                        *pointer = aligned_alloc(64, field.size_data);
#endif // _WIN32
                    }else{ // normal malloc
                        *pointer = realloc(*pointer, field.size_data);
                    }
                    fread(*pointer, field.size_data,1,inf);
                    *pointer_N = (size_t)field.size_data/fd.element_size;
                }
                goto next_field;
                break;
            case REB_STRING:
                // char** pointer = (char**)(base_address + fd.offset);
                // *pointer = realloc(*pointer, field.size_data); // TODO: memory needs to be freed somewhere else
                {
                    char* string = malloc(field.size_data);
                    fread(string, field.size_data,1,inf);
                    // HACK
                    reb_simulation_set_integrator(r, string);
                    free(string);
                    // /HACK
                }
                goto next_field;
                break;
            case REB_CHARP_LIST:
                {
                    size_t serialized_size = field.size_data;
                    char* serialized_strings = malloc(serialized_size);
                    fread(serialized_strings, serialized_size,1,inf);
                    // Process strings back into a list
                    size_t* pointer_N = (size_t*)fd.offset_N;
                    size_t current_pos = 0;
                    while (current_pos < serialized_size){
                        char* current_string = serialized_strings + current_pos;
                        // character count + NULL character + original pointer address
                        size_t current_string_len = strlen(current_string)+1+sizeof(char*);
                        current_pos += current_string_len;
                        // Add current_string to list
                        (*pointer_N)++;
                        *pointer = realloc(*pointer,sizeof(char*)*(*pointer_N));
                        (*(char***)pointer)[(*pointer_N)-1] = malloc(sizeof(char)*(current_string_len));
                        memcpy((*(char***)pointer)[(*pointer_N)-1], current_string, current_string_len);
                    }
                    free(serialized_strings);
                }
                goto next_field;
                break;
            case REB_OTHER:
                reb_simulation_warning(r, "Did not expect REB_OTHER here.");
                goto next_field;
                break;
            case REB_FIELD_NOT_FOUND:
                {
                    // We should never get here. If so, it's an unknown field id.
                    *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
                    int err = fseek(inf, field.size_data, SEEK_CUR);
                    if (err){
                        // Even worse, can't seek to end of field.
                        *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
                    }
                }
                goto finish_fields;
                break;
        }
        // If we're here then it was not a simple or pointer datatype. 
        // Can skip the iteration trough the descriptor list.
        // TODO:
        printf("shouldn't be hereeeee\n");
        break;
    }


finish_fields:
    // Some final initializations

    // Update pointers to simulation
    for (size_t l=0;l<r->N_var_config;l++){
        r->var_config[l].sim = r;
    }
    r->N_allocated = r->N; // This used to be different. Now only saving N.
    for (size_t l=0;l<r->N_allocated;l++){
        r->particles[l].ap = NULL;
        r->particles[l].sim = r;
#ifndef MPI
        // Restore names
        if (r->particles[l].name){
            int name_found = 0;
            for (size_t n=0;n<r->N_name_list;n++){
                char* original_pointer = *((char**)(r->name_list[n]+strlen(r->name_list[n])+1));
                if (r->particles[l].name==original_pointer){
                    name_found = 1;
                    r->particles[l].name = r->name_list[n];
                }
            }
            if (!name_found){
                reb_simulation_warning(r,"A name for a particle was not stored in the Simulationarchive.");
                r->particles[l].name = NULL;
            }
        }
#endif // MPI
    }
    for (size_t l=0;l<r->N_var;l++){
        r->particles_var[l].ap = NULL;
        r->particles_var[l].sim = r;
    }
}

struct reb_simulation* reb_binarydata_process_warnings(struct reb_simulation* r, enum REB_BINARYDATA_ERROR_CODE warnings){
    if (warnings & REB_BINARYDATA_ERROR_NOFILE){
        reb_simulation_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_WARNING_VERSION){
        reb_simulation_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if (warnings & REB_BINARYDATA_WARNING_POINTERS){
        reb_simulation_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_BINARYDATA_WARNING_PARTICLES){
        reb_simulation_warning(r,"Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if (warnings & REB_BINARYDATA_ERROR_FILENOTOPEN){
        reb_simulation_error(r,"Error while reading binary file (file was not open).");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_ERROR_OUTOFRANGE){
        reb_simulation_error(r,"Index out of range.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_ERROR_SEEK){
        reb_simulation_error(r,"Error while trying to seek file.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_WARNING_FIELD_UNKNOWN){
        reb_simulation_warning(r,"Unknown field found in binary file.");
    }
    if (warnings & REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR){
        reb_simulation_warning(r,"Custom integrator encountered in Simulationarchive. Call reb_simulation_set_integrator after the simulation is loaded to reset function pointers and initialize data.");
    }
    if (warnings & REB_BINARYDATA_ERROR_NOFILE){
        reb_simulation_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_ERROR_OLD){
        reb_simulation_error(r,"Reading old Simulationarchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_BINARYDATA_WARNING_CORRUPTFILE){
        reb_simulation_warning(r,"The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it.");
    }
    return r;
}

