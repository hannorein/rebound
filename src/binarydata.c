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
#include <stddef.h>
#include <inttypes.h>
#include "rebound.h"
#include "rebound_internal.h"
#include "particle.h"
#include "tools.h"
#include "tree.h"
#include "output.h"
#include "binarydata.h"
#include "simulationarchive.h"
#include "integrator_whfast512.h"
#include "integrator_janus.h"

// List of REBOUND parameters to be written to a file.
// Modify this list if you wish to input/output additional fields in the reb_simulation structure.
const struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_list[]= {
    { 0,  REB_DOUBLE,       "t",                            offsetof(struct reb_simulation, t), 0, 0},
    { 1,  REB_DOUBLE,       "G",                            offsetof(struct reb_simulation, G), 0, 0},
    { 2,  REB_DOUBLE,       "softening",                    offsetof(struct reb_simulation, softening), 0, 0},
    { 3,  REB_DOUBLE,       "dt",                           offsetof(struct reb_simulation, dt), 0, 0},
    { 4,  REB_SIZE_T,       "N",                            offsetof(struct reb_simulation, N), 0, 0},
    { 5,  REB_SIZE_T,       "N_var",                        offsetof(struct reb_simulation, N_var), 0, 0},
    // TODO: Add particles_var
    // 6 Used to be varconfig
    { 7,  REB_SIZE_T,       "N_active",                     offsetof(struct reb_simulation, N_active), 0, 0},
    { 8,  REB_INT,          "testparticle_type",            offsetof(struct reb_simulation, testparticle_type), 0, 0},
    { 10, REB_DOUBLE,       "opening_angle2",               offsetof(struct reb_simulation, opening_angle2), 0, 0},
    { 11, REB_INT,          "status",                       offsetof(struct reb_simulation, status), 0, 0},
    { 12, REB_INT,          "exact_finish_time",            offsetof(struct reb_simulation, exact_finish_time), 0, 0},
    { 13, REB_UINT,         "force_is_velocity_dependent",  offsetof(struct reb_simulation, force_is_velocity_dependent), 0, 0},
    { 14, REB_UINT,         "gravity_ignore_terms",         offsetof(struct reb_simulation, gravity_ignore_terms), 0, 0},
    { 15, REB_DOUBLE,       "output_timing_last",           offsetof(struct reb_simulation, output_timing_last), 0, 0},
    { 16, REB_INT,          "save_messages",                offsetof(struct reb_simulation, save_messages), 0, 0},
    { 17, REB_DOUBLE,       "exit_max_distance",            offsetof(struct reb_simulation, exit_max_distance), 0, 0},
    { 18, REB_DOUBLE,       "exit_min_distance",            offsetof(struct reb_simulation, exit_min_distance), 0, 0},
    { 19, REB_DOUBLE,       "usleep",                       offsetof(struct reb_simulation, usleep), 0, 0},
    { 20, REB_INT,          "track_energy_offset",          offsetof(struct reb_simulation, track_energy_offset), 0, 0},
    { 21, REB_DOUBLE,       "energy_offset",                offsetof(struct reb_simulation, energy_offset), 0, 0},
    // 22 was boxsize
    // 23 was boxsize_max
    { 24, REB_DOUBLE,       "root_size",                    offsetof(struct reb_simulation, root_size), 0, 0},
    // 25 was N_root
    { 26, REB_SIZE_T,       "N_root_x",                     offsetof(struct reb_simulation, N_root_x), 0, 0},
    { 27, REB_SIZE_T,       "N_root_y",                     offsetof(struct reb_simulation, N_root_y), 0, 0},
    { 28, REB_SIZE_T,       "N_root_z",                     offsetof(struct reb_simulation, N_root_z), 0, 0},
    { 29, REB_INT,          "N_ghost_x",                    offsetof(struct reb_simulation, N_ghost_x), 0, 0},
    { 30, REB_INT,          "N_ghost_y",                    offsetof(struct reb_simulation, N_ghost_y), 0, 0},
    { 31, REB_INT,          "N_ghost_z",                    offsetof(struct reb_simulation, N_ghost_z), 0, 0},
    // 32 was collision_resolve_keep_sorted
    { 33, REB_DOUBLE,       "minimum_collision_velocity",   offsetof(struct reb_simulation, minimum_collision_velocity), 0, 0},
    { 34, REB_DOUBLE,       "collisions_plog",              offsetof(struct reb_simulation, collisions_plog), 0, 0},
    { 36, REB_INT64,        "collisions_log_n",             offsetof(struct reb_simulation, collisions_log_n), 0, 0},
    { 37, REB_INT,          "calculate_megno",              offsetof(struct reb_simulation, calculate_megno), 0, 0},
    { 38, REB_DOUBLE,       "megno_Ys",                     offsetof(struct reb_simulation, megno_Ys), 0, 0},
    { 39, REB_DOUBLE,       "megno_Yss",                    offsetof(struct reb_simulation, megno_Yss), 0, 0},
    { 40, REB_DOUBLE,       "megno_cov_Yt",                 offsetof(struct reb_simulation, megno_cov_Yt), 0, 0},
    { 41, REB_DOUBLE,       "megno_var_t",                  offsetof(struct reb_simulation, megno_var_t), 0, 0},
    { 42, REB_DOUBLE,       "megno_mean_t",                 offsetof(struct reb_simulation, megno_mean_t), 0, 0},
    { 43, REB_DOUBLE,       "megno_mean_Y",                 offsetof(struct reb_simulation, megno_mean_Y), 0, 0},
    { 49, REB_DOUBLE,       "megno_initial_t",              offsetof(struct reb_simulation, megno_initial_t), 0, 0},
    { 44, REB_INT64,        "megno_n",                      offsetof(struct reb_simulation, megno_n), 0, 0},
    { 47, REB_DOUBLE,       "simulationarchive_auto_interval", offsetof(struct reb_simulation, simulationarchive_auto_interval), 0, 0},
    { 102, REB_DOUBLE,      "simulationarchive_auto_walltime", offsetof(struct reb_simulation, simulationarchive_auto_walltime), 0, 0},
    { 48, REB_DOUBLE,       "simulationarchive_next",       offsetof(struct reb_simulation, simulationarchive_next), 0, 0},
    { 50, REB_INT,          "collision",                    offsetof(struct reb_simulation, collision), 0, 0},
    { 51, REB_INT,          "integrator",                   offsetof(struct reb_simulation, integrator), 0, 0}, // Note: first element in structure
    { 52, REB_INT,          "boundary",                     offsetof(struct reb_simulation, boundary), 0, 0},
    { 53, REB_INT,          "gravity",                      offsetof(struct reb_simulation, gravity), 0, 0},
    { 54, REB_DOUBLE,       "OMEGA",                        offsetof(struct reb_simulation, OMEGA), 0, 0},
    { 55, REB_DOUBLE,       "OMEGAZ",                       offsetof(struct reb_simulation, OMEGAZ), 0, 0},
//    { 56, REB_DOUBLE,       "ri_sei.lastdt",                offsetof(struct reb_simulation, ri_sei.lastdt), 0, 0},
//    { 57, REB_DOUBLE,       "ri_sei.sindt",                 offsetof(struct reb_simulation, ri_sei.sindt), 0, 0},
//    { 58, REB_DOUBLE,       "ri_sei.tandt",                 offsetof(struct reb_simulation, ri_sei.tandt), 0, 0},
//    { 59, REB_DOUBLE,       "ri_sei.sindtz",                offsetof(struct reb_simulation, ri_sei.sindtz), 0, 0},
//    { 60, REB_DOUBLE,       "ri_sei.tandtz",                offsetof(struct reb_simulation, ri_sei.tandtz), 0, 0},
    { 69, REB_DOUBLE,       "ri_ias15.epsilon",             offsetof(struct reb_simulation, ri_ias15.epsilon), 0, 0},
    { 70, REB_DOUBLE,       "ri_ias15.min_dt",              offsetof(struct reb_simulation, ri_ias15.min_dt), 0, 0},
    { 71, REB_UINT,         "ri_ias15.adaptive_mode",       offsetof(struct reb_simulation, ri_ias15.adaptive_mode), 0, 0},
    { 72, REB_UINT64,        "ri_ias15.iterations_max_exceeded", offsetof(struct reb_simulation, ri_ias15.iterations_max_exceeded), 0, 0},
    { 85, REB_POINTER,      "particles",                    offsetof(struct reb_simulation, particles), offsetof(struct reb_simulation, N), sizeof(struct reb_particle)},
    { 403, REB_POINTER,     "particles_var",                offsetof(struct reb_simulation, particles_var), offsetof(struct reb_simulation, N_var), sizeof(struct reb_particle)},
    { 86, REB_POINTER,      "var_config",                   offsetof(struct reb_simulation, var_config), offsetof(struct reb_simulation, N_var_config), sizeof(struct reb_variational_configuration)},
    { 87, REB_OTHER,        "functionpointers", 0, 0, 0},
    { 89, REB_POINTER,      "ri_ias15.at",                  offsetof(struct reb_simulation, ri_ias15.at), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 90, REB_POINTER,      "ri_ias15.x0",                  offsetof(struct reb_simulation, ri_ias15.x0), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 91, REB_POINTER,      "ri_ias15.v0",                  offsetof(struct reb_simulation, ri_ias15.v0), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 92, REB_POINTER,      "ri_ias15.a0",                  offsetof(struct reb_simulation, ri_ias15.a0), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 93, REB_POINTER,      "ri_ias15.csx",                 offsetof(struct reb_simulation, ri_ias15.csx), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 94, REB_POINTER,      "ri_ias15.csv",                 offsetof(struct reb_simulation, ri_ias15.csv), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 95, REB_POINTER,      "ri_ias15.csa0",                offsetof(struct reb_simulation, ri_ias15.csa0), offsetof(struct reb_simulation, ri_ias15.N_allocated), sizeof(double)},
    { 96, REB_DP7,          "ri_ias15.g",                   offsetof(struct reb_simulation, ri_ias15.g), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    { 97, REB_DP7,          "ri_ias15.b",                   offsetof(struct reb_simulation, ri_ias15.b), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    { 98, REB_DP7,          "ri_ias15.csb",                 offsetof(struct reb_simulation, ri_ias15.csb), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    { 99, REB_DP7,          "ri_ias15.e",                   offsetof(struct reb_simulation, ri_ias15.e), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    { 100, REB_DP7,         "ri_ias15.br",                  offsetof(struct reb_simulation, ri_ias15.br), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    { 101, REB_DP7,         "ri_ias15.er",                  offsetof(struct reb_simulation, ri_ias15.er), offsetof(struct reb_simulation, ri_ias15.N_allocated), 7*sizeof(double)},
    //{ 107, REB_INT,         "visualization",                offsetof(struct reb_simulation, visualization), 0, 0},
    { 112, REB_POINTER,     "ri_janus.p_int",               offsetof(struct reb_simulation, ri_janus.p_int), offsetof(struct reb_simulation, ri_janus.N_allocated), sizeof(struct reb_particle_int)},
    { 113, REB_DOUBLE,      "ri_janus.scale_pos",           offsetof(struct reb_simulation, ri_janus.scale_pos), 0, 0},
    { 114, REB_DOUBLE,      "ri_janus.scale_vel",           offsetof(struct reb_simulation, ri_janus.scale_vel), 0, 0},
    { 115, REB_UINT,        "ri_janus.order",               offsetof(struct reb_simulation, ri_janus.order), 0, 0},
    { 116, REB_UINT,        "ri_janus.recalculate_integer_coordinates_this_timestep", offsetof(struct reb_simulation, ri_janus.recalculate_integer_coordinates_this_timestep), 0, 0},
    { 118, REB_DOUBLE,      "ri_mercurius.r_crit_hill",     offsetof(struct reb_simulation, ri_mercurius.r_crit_hill), 0, 0},
    { 119, REB_UINT,        "ri_mercurius.safe_mode",       offsetof(struct reb_simulation, ri_mercurius.safe_mode), 0, 0},
    { 120, REB_UINT,        "ri_mercurius.is_synchronized", offsetof(struct reb_simulation, ri_mercurius.is_synchronized), 0, 0},
    { 122, REB_POINTER,     "ri_mercurius.dcrit",           offsetof(struct reb_simulation, ri_mercurius.dcrit), offsetof(struct reb_simulation, ri_mercurius.N_allocated_dcrit), sizeof(double)},
    { 123, REB_UINT,        "ri_mercurius.recalculate_coordinates_this_timestep", offsetof(struct reb_simulation, ri_mercurius.recalculate_coordinates_this_timestep), 0, 0},
    { 125, REB_INT,         "simulationarchive_version",    offsetof(struct reb_simulation, simulationarchive_version), 0, 0},
    { 126, REB_DOUBLE,      "walltime",                     offsetof(struct reb_simulation, walltime), 0, 0},
    { 127, REB_DOUBLE,      "walltime_last_steps",          offsetof(struct reb_simulation, walltime_last_steps), 0, 0},
    { 130, REB_UINT32,      "python_unit_l",                offsetof(struct reb_simulation, python_unit_l), 0, 0},
    { 131, REB_UINT32,      "python_unit_m",                offsetof(struct reb_simulation, python_unit_m), 0, 0},
    { 132, REB_UINT32,      "python_unit_t",                offsetof(struct reb_simulation, python_unit_t), 0, 0},
    { 133, REB_VEC3D,       "ri_mercurius.com_pos",         offsetof(struct reb_simulation, ri_mercurius.com_pos), 0, 0},
    { 134, REB_VEC3D,       "ri_mercurius.com_vel",         offsetof(struct reb_simulation, ri_mercurius.com_vel), 0, 0},
    { 135, REB_UINT64,      "simulationarchive_auto_step",  offsetof(struct reb_simulation, simulationarchive_auto_step), 0, 0},
    { 136, REB_UINT64,      "simulationarchive_next_step",  offsetof(struct reb_simulation, simulationarchive_next_step), 0, 0},
    { 137, REB_UINT64,      "steps_done",                   offsetof(struct reb_simulation, steps_done), 0, 0},
    { 140, REB_UINT,        "ri_saba.safe_mode",            offsetof(struct reb_simulation, ri_saba.safe_mode), 0, 0},
    { 141, REB_UINT,        "ri_saba.is_synchronized",      offsetof(struct reb_simulation, ri_saba.is_synchronized), 0, 0},
    { 145, REB_DOUBLE,      "dt_last_done",                 offsetof(struct reb_simulation, dt_last_done), 0, 0},
    { 146, REB_INT,         "ri_saba.type",                 offsetof(struct reb_simulation, ri_saba.type), 0, 0},
    { 147, REB_UINT,        "ri_saba.keep_unsynchronized",  offsetof(struct reb_simulation, ri_saba.keep_unsynchronized), 0, 0},
    { 148, REB_INT,         "ri_eos.phi0",                  offsetof(struct reb_simulation, ri_eos.phi0), 0, 0},
    { 149, REB_INT,         "ri_eos.phi1",                  offsetof(struct reb_simulation, ri_eos.phi1), 0, 0},
    { 150, REB_UINT,        "ri_eos.n",                     offsetof(struct reb_simulation, ri_eos.n), 0, 0},
    { 151, REB_UINT,        "ri_eos.safe_mode",             offsetof(struct reb_simulation, ri_eos.safe_mode), 0, 0},
    { 152, REB_UINT,        "ri_eos.is_synchronized",       offsetof(struct reb_simulation, ri_eos.is_synchronized), 0, 0},
    { 154, REB_UINT,        "rand_seed",                    offsetof(struct reb_simulation, rand_seed), 0, 0},
    { 155, REB_INT,         "testparticle_hidewarnings",    offsetof(struct reb_simulation, testparticle_hidewarnings), 0, 0},
    { 156, REB_DOUBLE,      "ri_bs.eps_abs",                offsetof(struct reb_simulation, ri_bs.eps_abs), 0, 0},
    { 157, REB_DOUBLE,      "ri_bs.eps_rel",                offsetof(struct reb_simulation, ri_bs.eps_rel), 0, 0},
    { 158, REB_DOUBLE,      "ri_bs.min_dt",                 offsetof(struct reb_simulation, ri_bs.min_dt), 0, 0},
    { 159, REB_DOUBLE,      "ri_bs.max_dt",                 offsetof(struct reb_simulation, ri_bs.max_dt), 0, 0},
    { 160, REB_INT,         "ri_bs.first_or_last_step",     offsetof(struct reb_simulation, ri_bs.first_or_last_step), 0, 0},
    { 161, REB_INT,         "ri_bs.previous_rejected",      offsetof(struct reb_simulation, ri_bs.previous_rejected), 0, 0},
    { 162, REB_INT,         "ri_bs.target_iter",            offsetof(struct reb_simulation, ri_bs.target_iter), 0, 0},
    { 164, REB_POINTER_FIXED_SIZE, "display_settings",      offsetof(struct reb_simulation, display_settings), 0, sizeof(struct reb_display_settings)},
    { 165, REB_DOUBLE,      "ri_trace.r_crit_hill",         offsetof(struct reb_simulation, ri_trace.r_crit_hill), 0, 0},
    // TRACE Pericenter conditions used to have ids 166 - 168. Do not reuse.
    { 169, REB_DOUBLE,      "ri_trace.peri_crit_eta",       offsetof(struct reb_simulation, ri_trace.peri_crit_eta), 0, 0},
    { 170, REB_INT, "ri_trace.peri_mode", offsetof(struct reb_simulation, ri_trace.peri_mode), 0, 0},
    //    { 163, REB_INT,         "var_rescale_warning", offsetof(struct reb_simulation, var_rescale_warning), 0, 0},
    // TES Variables used to have ids 300 - 388. Do not reuse. 
    { 390, REB_UINT,        "ri_whfast512.keep_unsynchronized", offsetof(struct reb_simulation, ri_whfast512.keep_unsynchronized), 0, 0},
    { 391, REB_UINT,        "ri_whfast512.is_synchronized", offsetof(struct reb_simulation, ri_whfast512.is_synchronized), 0, 0},
    { 392, REB_UINT,        "ri_whfast512.gr_potential",    offsetof(struct reb_simulation, ri_whfast512.gr_potential), 0, 0},
    { 394, REB_POINTER_ALIGNED, "ri_whfast512.pjh",         offsetof(struct reb_simulation, ri_whfast512.p_jh), offsetof(struct reb_simulation, ri_whfast512.N_allocated), sizeof(struct reb_particle_avx512)},
    // 396, 397 used to be max_radius0 and max_radius1
    { 398, REB_UINT,        "ri_whfast512.N_systems",       offsetof(struct reb_simulation, ri_whfast512.N_systems), 0, 0},
    { 399, REB_PARTICLE4,   "ri_whfast512.pjh0",            offsetof(struct reb_simulation, ri_whfast512.p_jh0), 0, 0},
    { 402, REB_CHARP_LIST,  "name_list",                    offsetof(struct reb_simulation, name_list), offsetof(struct reb_simulation, N_name_list), 0},
    // 403  particles_var
    { 1329743186, REB_OTHER,"header", 0, 0, 0},
    { 9998, REB_OTHER,      "sablob", 0, 0, 0},
    { 9999, REB_FIELD_END,  "end", 0, 0, 0}
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

struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_type(uint32_t type){
    int i=-1;
    do{
        i++;
        if (reb_binarydata_field_descriptor_list[i].type==type){
            return reb_binarydata_field_descriptor_list[i];
        }
    } while (reb_binarydata_field_descriptor_list[i].dtype!=REB_FIELD_END);
    struct reb_binarydata_field_descriptor bfd = {0};
    bfd.dtype = REB_FIELD_NOT_FOUND;
    return bfd;
}

struct reb_binarydata_field_descriptor reb_binarydata_field_descriptor_for_name(const char* name){
    int i=-1;
    do{
        i++;
        if (strcmp(reb_binarydata_field_descriptor_list[i].name, name)==0){
            return reb_binarydata_field_descriptor_list[i];
        }
    } while (reb_binarydata_field_descriptor_list[i].dtype!=REB_FIELD_END);
    struct reb_binarydata_field_descriptor bfd = {0};
    bfd.dtype = REB_FIELD_NOT_FOUND;
    return bfd;
}

// Helper function to print out binary data in human readable form.
static void asprintf_reb_type(char** buf, enum REB_BINARYDATA_DTYPE dtype, char* pointer, size_t dsize){
    char* newbuf = NULL;
    switch (dtype){
        case REB_DOUBLE:
            asprintf(&newbuf,"%e",*(double*)(pointer));
            break;
        case REB_INT:
            asprintf(&newbuf,"%d",*(int*)(pointer));
            break;
        case REB_SIZE_T:
            asprintf(&newbuf,"%zu",*(size_t*)(pointer));
            break;
        case REB_UINT:
            asprintf(&newbuf,"%u",*(unsigned int*)(pointer));
            break;
        case REB_UINT32:
            asprintf(&newbuf,"%" PRIu32,*(uint32_t*)(pointer)); // PRIu32 defined in inttypes.h
            break;
        case REB_INT64:
            asprintf(&newbuf,"%" PRId64,*(int64_t*)(pointer));
            break;
        case REB_UINT64:
            asprintf(&newbuf,"%" PRIu64,*(uint64_t*)(pointer));
            break;
        default:
            asprintf(&newbuf,"(%zu bytes, values not printed)", dsize);
            break;
    }
    if (buf){
        *buf = realloc(*buf, strlen(*buf) + strlen(newbuf) + sizeof(char));
        strcat(*buf,newbuf);
    }else{
        printf("%s",newbuf);
    }
    free(newbuf);
}

// Compares two simulations in buffers.
// Returns 0 if the buffers contain the same simulation data. 
// Supports different output options.
int reb_binarydata_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option){
    if (!buf1 || !buf2 || size1<64 || size2<64){
        printf("Cannot read input buffers.\n");
        return 0;
    }

    int are_different = 0;

    if (output_option==0){
        *bufp = NULL;
        *sizep = 0;
    }
    if (output_option==3){
        *bufp = malloc(sizeof(char));
        *bufp[0] = '\0';
    }
    size_t allocatedsize = 0;

    // Header.
    if(memcmp(buf1,buf2,64)!=0 && output_option==1){
        printf("Header in binary files are different.\n");
    }

    size_t pos1 = 64;
    size_t pos2 = 64;

    struct reb_binarydata_field_descriptor fd_end = reb_binarydata_field_descriptor_for_name("end");

    while(1){
        if (pos1+sizeof(struct reb_binarydata_field)>size1) break;
        struct reb_binarydata_field field1;
        memcpy(&field1, buf1+pos1, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos1 += sizeof(struct reb_binarydata_field);
        if (field1.type==fd_end.type){
            break;
        }
        if (pos2+sizeof(struct reb_binarydata_field)>size2) pos2 = 64;
        struct reb_binarydata_field field2;
        memcpy(&field2, buf2+pos2, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos2 += sizeof(struct reb_binarydata_field);

        // Fields might not be in the same order.
        if (field1.type!=field2.type){
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
                if(field2.type==fd_end.type){
                    notfound = 1;
                    break;
                }
                if (field2.type==field1.type){
                    break; // found!!
                }else{
                    pos2 += field2.size; //skip
                }
            };
            if (notfound == 1){
                pos1 += field1.size; // For next search
                pos2 = 64;           // For next search
                field1.size = 0;     // Output field with size 0
                are_different = 1.;
                if (output_option==0){
                    write_to_stream(bufp, &allocatedsize, sizep, &field1,sizeof(struct reb_binarydata_field));
                }else if (output_option==1 || output_option==3){
                    const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_type(field1.type);
                    char* buf;
#ifndef _WIN32
                    asprintf(&buf, "%s:\n\033[31m< ",fd.name);
#else // _WIN32
                    asprintf(&buf, "%s:\n< ",fd.name);
#endif // _WIN32
                    if (bufp){
                        *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                        strcat(*bufp,buf);
                    }else{
                        printf("%s",buf);
                    }
                    free(buf);
                    asprintf_reb_type(bufp, fd.dtype, buf1+pos1, field1.size);
#ifndef _WIN32
                    asprintf(&buf, "\033[0m\n");
#else // _WIN32
                    asprintf(&buf, "\n");
#endif // _WIN32
                    if (bufp){
                        *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                        strcat(*bufp,buf);
                    }else{
                        printf("%s",buf);
                    }
                    free(buf);
                }
                field1.size = 0;
                continue;
            }
        }
        // Can assume field1.type == field2.type from here on
        if (pos1+field1.size>size1) printf("Corrupt binary file buf1.\n");
        if (pos2+field2.size>size2) printf("Corrupt binary file buf2.\n");
        int fields_differ = 0;
        if (field1.size==field2.size){
            if (strcmp(reb_binarydata_field_descriptor_for_type(field1.type).name, "particles")==0){
                struct reb_particle* pb1 = (struct reb_particle*)(buf1+pos1);
                struct reb_particle* pb2 = (struct reb_particle*)(buf2+pos2);
                for (size_t i=0;i<field1.size/sizeof(struct reb_particle);i++){
                    struct reb_particle p1;
                    struct reb_particle p2;
                    memcpy(&p1, pb1+i, sizeof(struct reb_particle)); // need copy because of 8 byte alignment requirement
                    memcpy(&p2, pb2+i, sizeof(struct reb_particle)); // need copy because of 8 byte alignment requirement
                    fields_differ |= reb_particle_cmp(p1,p2);
                }
            }else{
                if (memcmp(buf1+pos1,buf2+pos2,field1.size)!=0){
                    fields_differ = 1;
                }
            }
        }else{
            fields_differ = 1;
        }
        if(fields_differ){
            if (strncmp(reb_binarydata_field_descriptor_for_type(field1.type).name, "walltime",8)!=0){
                // Ignore the walltime fields, but only for the return value (print it out)
                // Typically we do not care about this field when comparing simulations.
                are_different = 1.;
            }
            if (output_option==0){
                write_to_stream(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, &allocatedsize, sizep, buf2+pos2,field2.size);
            }else if (output_option==1 || output_option==3){
                const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_type(field1.type);
                char* buf;
#ifndef _WIN32
                asprintf(&buf, "%s:\n\033[31m< ",fd.name);
#else // _WIN32
                asprintf(&buf, "%s:\n< ",fd.name);
#endif // _WIN32
                if (bufp){
                    *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                    strcat(*bufp,buf);
                }else{
                    printf("%s",buf);
                }
                free(buf);
                asprintf_reb_type(bufp, fd.dtype, buf1+pos1, field1.size);
#ifndef _WIN32
                asprintf(&buf, "\033[0m\n---\n\033[32m> ");
#else // _WIN32
                asprintf(&buf, "\n---\n> ");
#endif // _WIN32
                if (bufp){
                    *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                    strcat(*bufp,buf);
                }else{
                    printf("%s",buf);
                }
                free(buf);
                asprintf_reb_type(bufp, fd.dtype, buf2+pos2, field2.size);
#ifndef _WIN32
                asprintf(&buf, "\033[0m\n");
#else // _WIN32
                asprintf(&buf, "\n");
#endif // _WIN32
                if (bufp){
                    *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                    strcat(*bufp,buf);
                }else{
                    printf("%s",buf);
                }
                free(buf);
            }
        }
        pos1 += field1.size;
        pos2 += field2.size;
    }
    // Search for fields which are present in buf2 but not in buf1
    pos1 = 64;
    pos2 = 64;
    while(1){
        if (pos2+sizeof(struct reb_binarydata_field)>size2) break;
        struct reb_binarydata_field field2;
        memcpy(&field2, buf2+pos2, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos2 += sizeof(struct reb_binarydata_field);
        if (field2.type==fd_end.type){
            break;
        }
        if (pos1+sizeof(struct reb_binarydata_field)>size1) pos1 = 64;
        struct reb_binarydata_field field1;
        memcpy(&field1, buf1+pos1, sizeof(struct reb_binarydata_field)); // need copy because of 8 byte alignment requirement
        pos1 += sizeof(struct reb_binarydata_field);

        if (field1.type==field2.type){
            // Not a new field. Skip.
            pos1 += field1.size;
            pos2 += field2.size;
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
            if(field1.type==fd_end.type){
                notfound = 1;
                break;
            }
            if (field2.type==field1.type){
                break; // found it, not new
            }else{
                // not found, try next
                pos1 += field1.size;
            }
        };
        if (notfound == 0){
            // Not a new field. Skip.
            pos1 = 64;
            pos2 += field2.size;
            continue;
        }

        are_different = 1.;
        if (output_option==0){
            write_to_stream(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binarydata_field));
            write_to_stream(bufp, &allocatedsize, sizep, buf2+pos2,field2.size);
        }else if (output_option==1 || output_option==3){
            const struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_for_type(field2.type);
            char* buf;
#ifndef _WIN32
            asprintf(&buf, "%s:\n\033[32m> ",fd.name);
#else // _WIN32
            asprintf(&buf, "%s:\n> ",fd.name);
#endif // _WIN32
            if (bufp){
                *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                strcat(*bufp,buf);
            }else{
                printf("%s",buf);
            }
            asprintf_reb_type(bufp, fd.dtype, buf2+pos2, field2.size);
#ifndef _WIN32
            asprintf(&buf, "\033[0m\n");
#else // _WIN32
            asprintf(&buf, "\n");
#endif // _WIN32
            if (bufp){
                *bufp = realloc(*bufp, strlen(*bufp) + strlen(buf) + sizeof(char));
                strcat(*bufp,buf);
            }else{
                printf("%s",buf);
            }
        }
        pos1 = 64;
        pos2 += field2.size;
    }

    return are_different;
}


// Macro to write a single field to a binary file.
// Memset forces padding to be set to 0 (not necessary but
// helps when comparing binary files)
#define WRITE_FIELD_TYPE(typen, value, length) {\
    struct reb_binarydata_field field;\
    memset(&field,0,sizeof(struct reb_binarydata_field));\
    field.type = typen;\
    field.size = (length);\
    write_to_stream(bufp, &allocatedsize, sizep, &field,sizeof(struct reb_binarydata_field));\
    write_to_stream(bufp, &allocatedsize, sizep, value,field.size);\
}

// Serializes a simulation to a buffer
void reb_binarydata_simulation_to_stream(struct reb_simulation* r, char** bufp, size_t* sizep){
    if (r->simulationarchive_version<3){
        reb_simulation_error(r, "Simulationarchives with version < 3 are no longer supported.\n");
    }
    size_t allocatedsize = 0;
    *bufp = NULL;
    *sizep = 0;

    // Output header.
    char header[64] = "\0";
    int cwritten = sprintf(header,"REBOUND Binary File. Version: %s",reb_version_str);
    snprintf(header+cwritten+1,64-cwritten-1,"%s",reb_githash_str);
    write_to_stream(bufp, &allocatedsize, sizep, header,sizeof(char)*64);

    // Compress data if possible
    // This does not affect future calculation, but might trigger a realloc.
    if (r->ri_ias15.N_allocated > 3*r->N){
        r->ri_ias15.N_allocated = 3*r->N;
    }
    /// Output all fields
    int i=0;
    while (reb_binarydata_field_descriptor_list[i].dtype!=REB_FIELD_END){
        int dtype = reb_binarydata_field_descriptor_list[i].dtype;
        // Simple data types:
        if (dtype == REB_DOUBLE || dtype == REB_INT || dtype == REB_UINT || dtype == REB_UINT32
                || dtype == REB_INT64 || dtype == REB_UINT64 || dtype == REB_PARTICLE 
                || dtype == REB_PARTICLE4 || dtype == REB_VEC3D || dtype == REB_SIZE_T){
            struct reb_binarydata_field field;
            memset(&field,0,sizeof(struct reb_binarydata_field));
            field.type = reb_binarydata_field_descriptor_list[i].type;
            switch (dtype){
                case REB_DOUBLE: 
                    field.size = sizeof(double);
                    break;
                case REB_INT: 
                    field.size = sizeof(int);
                    break;
                case REB_SIZE_T: 
                    field.size = sizeof(size_t);
                    break;
                case REB_UINT: 
                    field.size = sizeof(unsigned int);
                    break;
                case REB_UINT32: 
                    field.size = sizeof(uint32_t);
                    break;
                case REB_INT64:
                    field.size = sizeof(int64_t);
                    break;
                case REB_UINT64:
                    field.size = sizeof(uint64_t);
                    break;
                case REB_VEC3D:
                    field.size = sizeof(struct reb_vec3d);
                    break;
                case REB_PARTICLE:
                    field.size = sizeof(struct reb_particle);
                    break;
                case REB_PARTICLE4:
                    field.size = 4*sizeof(struct reb_particle);
                    break;
            }
            write_to_stream(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binarydata_field));
            char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
            write_to_stream(bufp, &allocatedsize, sizep, pointer, field.size);
        }
        // Pointer data types
        if (dtype == REB_POINTER || dtype == REB_POINTER_ALIGNED ){
            struct reb_binarydata_field field;
            memset(&field,0,sizeof(struct reb_binarydata_field));
            field.type = reb_binarydata_field_descriptor_list[i].type;
            size_t* pointer_N = (size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binarydata_field_descriptor_list[i].element_size;

            if (field.size){
                write_to_stream(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binarydata_field));
                char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                pointer = *(char**)pointer;
                write_to_stream(bufp, &allocatedsize, sizep, pointer, field.size);
            }
        }
        // Pointer with a fixed size
        if (dtype == REB_POINTER_FIXED_SIZE ){
            struct reb_binarydata_field field;
            memset(&field,0,sizeof(struct reb_binarydata_field));
            field.type = reb_binarydata_field_descriptor_list[i].type;
            field.size = reb_binarydata_field_descriptor_list[i].element_size;

            char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
            pointer = *(char**)pointer;
            if (pointer){
                write_to_stream(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binarydata_field));
                write_to_stream(bufp, &allocatedsize, sizep, pointer, field.size);
            }
        }
        if (dtype == REB_CHARP_LIST ){
            struct reb_binarydata_field field;
            memset(&field,0,sizeof(struct reb_binarydata_field));
            field.type = reb_binarydata_field_descriptor_list[i].type;
            size_t N_list = *((size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N));
            char*** list_p = (char***)((char*)r + reb_binarydata_field_descriptor_list[i].offset);
            size_t serialized_size = 0;
            for (size_t i=0; i<N_list; i++){
                // character count + NULL character + original pointer
                serialized_size += strlen((*list_p)[i])+1+sizeof(char*);
            }
            field.size = sizeof(char)*serialized_size;

            if (field.size){
                // This pointer arithmetic will fail on 32 bit architectures.
                write_to_stream(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binarydata_field));
                for (size_t i=0; i<N_list; i++){
                    write_to_stream(bufp, &allocatedsize, sizep, (*list_p)[i], strlen((*list_p)[i])+1);
                    write_to_stream(bufp, &allocatedsize, sizep, &((*list_p)[i]), sizeof(char*));
                }
            }
        }
        // Special datatype for IAS15. Similar to POINTER
        if (dtype == REB_DP7 ){
            struct reb_binarydata_field field;
            memset(&field,0,sizeof(struct reb_binarydata_field));
            field.type = reb_binarydata_field_descriptor_list[i].type;
            size_t* pointer_N = (size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N);
            field.size = (*pointer_N) * reb_binarydata_field_descriptor_list[i].element_size;

            if (field.size){
                write_to_stream(bufp, &allocatedsize, sizep, &field, sizeof(struct reb_binarydata_field));
                char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                struct reb_dp7* dp7 = (struct reb_dp7*)pointer;
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p0,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p1,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p2,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p3,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p4,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p5,field.size/7);
                write_to_stream(bufp, &allocatedsize, sizep, dp7->p6,field.size/7);
            }
        }
        i++;
    }

    // Write function pointer warning flag
    int functionpointersused = 0;
    if (r->coefficient_of_restitution ||
            r->collision_resolve ||
            r->additional_forces ||
            r->heartbeat ||
            r->ri_trace.S ||
            r->ri_trace.S_peri ||
            r->post_timestep_modifications ||
            r->free_particle_ap){
        functionpointersused = 1;
    }

    struct reb_binarydata_field_descriptor fd_fp = reb_binarydata_field_descriptor_for_name("functionpointers");
    struct reb_binarydata_field field_functionp;
    memset(&field_functionp,0,sizeof(struct reb_binarydata_field));
    field_functionp.type = fd_fp.type; 
    field_functionp.size = sizeof(int);
    write_to_stream(bufp, &allocatedsize, sizep, &field_functionp, sizeof(struct reb_binarydata_field));
    write_to_stream(bufp, &allocatedsize, sizep, &functionpointersused, field_functionp.size);

    int end_null = 0;

    struct reb_binarydata_field_descriptor fd_end = reb_binarydata_field_descriptor_for_name("end");
    WRITE_FIELD_TYPE(fd_end.type, &end_null, 0);
    struct reb_simulationarchive_blob blob = {0};
    write_to_stream(bufp, &allocatedsize, sizep, &blob, sizeof(struct reb_simulationarchive_blob));
}

// Read field data into simulation from file or memory buffer.
void reb_binarydata_input_fields(struct reb_simulation* r, FILE* inf, enum REB_BINARYDATA_ERROR_CODE* warnings){
    struct reb_binarydata_field field;
    // A few fields need special treatment. Find their descriptors first.
    struct reb_binarydata_field_descriptor fd_header = reb_binarydata_field_descriptor_for_name("header");
    struct reb_binarydata_field_descriptor fd_end = reb_binarydata_field_descriptor_for_name("end");
    struct reb_binarydata_field_descriptor fd_functionpointers = reb_binarydata_field_descriptor_for_name("functionpointers");

next_field:
    // Loop over all fields
    while(1){

        int numread = (int)fread(&field,sizeof(struct reb_binarydata_field),1,inf);
        if (numread<1){
            goto finish_fields; // End of file
        }
        if (field.type==fd_end.type){
            goto finish_fields; // End of snapshot
        }
        int i=0;

        // Loop over field descriptor list. Simple datatypes and pointers will be read in this loop.
        while (reb_binarydata_field_descriptor_list[i].dtype!=REB_FIELD_END){
            struct reb_binarydata_field_descriptor fd = reb_binarydata_field_descriptor_list[i];
            if (fd.type==field.type){
                // Read simple data types
                if (fd.dtype == REB_DOUBLE || fd.dtype == REB_INT || fd.dtype == REB_UINT 
                        || fd.dtype == REB_UINT32 || fd.dtype == REB_INT64 
                        || fd.dtype == REB_UINT64 || fd.dtype == REB_PARTICLE 
                        || fd.dtype == REB_PARTICLE4 || fd.dtype == REB_VEC3D || fd.dtype == REB_SIZE_T ){
                    char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                    fread(pointer, field.size, 1, inf);
                    goto next_field;
                }
                // Read a pointer data type. 
                // 1) reallocate memory
                // 2) read data into memory
                // 3) set N_allocated variable
                if (fd.dtype == REB_POINTER || fd.dtype == REB_POINTER_ALIGNED){
                    if (field.size % reb_binarydata_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field.");
                    }
                    char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                    if (fd.dtype == REB_POINTER_ALIGNED){
                        if (*(char**)pointer) free(*(char**)pointer);
#if defined(_WIN32) || !defined(AVX512)
                        // WHFast512 not supported on Windows!
                        *(char**)pointer = malloc(sizeof(struct reb_particle_avx512));
#else 
                        *(char**)pointer = aligned_alloc(64,sizeof(struct reb_particle_avx512));
#endif // _WIN32
                    }else{ // normal malloc
                        *(char**)pointer = realloc(*(char**)pointer, field.size);
                    }
                    fread(*(char**)pointer, field.size,1,inf);

                    size_t* pointer_N = (size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N);
                    *pointer_N = (size_t)field.size/reb_binarydata_field_descriptor_list[i].element_size;

                    goto next_field;
                }
                if (fd.dtype == REB_POINTER_FIXED_SIZE){
                    if (field.size != reb_binarydata_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field (fixed pointer size).");
                    }
                    char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                    *(char**)pointer = realloc(*(char**)pointer, field.size);
                    fread(*(char**)pointer, field.size,1,inf);

                    goto next_field;
                }
                if (fd.dtype == REB_CHARP_LIST){
                    size_t serialized_size = field.size;
                    char* serialized_strings = malloc(serialized_size);
                    fread(serialized_strings, serialized_size,1,inf);
                    // Process strings back into a list
                    char*** pointer = (char***)((char*)r + reb_binarydata_field_descriptor_list[i].offset);
                    size_t* pointer_N = (size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N);
                    size_t current_pos = 0;
                    while (current_pos < serialized_size){
                        char* current_string = serialized_strings + current_pos;
                        // character count + NULL character + original pointer address
                        size_t current_string_len = strlen(current_string)+1+sizeof(char*);
                        current_pos += current_string_len;
                        // Add current_string to list
                        (*pointer_N)++;
                        *pointer = realloc(*pointer,sizeof(char*)*(*pointer_N));
                        (*pointer)[(*pointer_N)-1] = malloc(sizeof(char)*(current_string_len));
                        memcpy((*pointer)[(*pointer_N)-1], current_string, current_string_len);
                    }
                    free(serialized_strings);
                    goto next_field;
                }
                // Special datatype for ias15. Similar to REB_POINTER. 
                if (fd.dtype == REB_DP7){
                    if (field.size % reb_binarydata_field_descriptor_list[i].element_size){
                        reb_simulation_warning(r, "Inconsistent size encountered in binary field.");
                    }
                    char* pointer = (char*)r + reb_binarydata_field_descriptor_list[i].offset;
                    struct reb_dp7* dp7 = (struct reb_dp7*)pointer;

                    dp7->p0 = realloc(dp7->p0,field.size/7);
                    dp7->p1 = realloc(dp7->p1,field.size/7);
                    dp7->p2 = realloc(dp7->p2,field.size/7);
                    dp7->p3 = realloc(dp7->p3,field.size/7);
                    dp7->p4 = realloc(dp7->p4,field.size/7);
                    dp7->p5 = realloc(dp7->p5,field.size/7);
                    dp7->p6 = realloc(dp7->p6,field.size/7);
                    fread(dp7->p0, field.size/7, 1, inf);
                    fread(dp7->p1, field.size/7, 1, inf);
                    fread(dp7->p2, field.size/7, 1, inf);
                    fread(dp7->p3, field.size/7, 1, inf);
                    fread(dp7->p4, field.size/7, 1, inf);
                    fread(dp7->p5, field.size/7, 1, inf);
                    fread(dp7->p6, field.size/7, 1, inf);

                    size_t* pointer_N = (size_t*)((char*)r + reb_binarydata_field_descriptor_list[i].offset_N);
                    *pointer_N = (size_t)field.size/reb_binarydata_field_descriptor_list[i].element_size;

                    goto next_field;
                }
                // If we're here then it was not a simple or pointer datatype. 
                // Can skip the iteration trough the descriptor list.
                break;
            }
            i++;
        }

        // Fields with types that require special handling
        if (field.type == fd_functionpointers.type){
            // Warning for when function pointers were used. 
            // No effect on simulation.
            int fpwarn;
            fread(&fpwarn, field.size,1,inf);
            if (fpwarn && warnings){
                *warnings |= REB_BINARYDATA_WARNING_POINTERS;
            }
            goto next_field;
        }
        if (field.type == fd_header.type){
            // Check header.
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

        // We should never get here. If so, it's an unknown field type.
        *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
        int err = fseek(inf, field.size, SEEK_CUR);
        if (err){
            // Even worse, can't seek to end of field.
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
        }
    } 

finish_fields:
    // Some final initializations

    // Find integrator
    {
        int integrator_found = 0;
        for (size_t i = 0; i<reb_integrators_available_N;i++){
            if (reb_integrators_available[i]->id == r->integrator.id){
                r->integrator = *reb_integrators_available[i];
                integrator_found = 1;
                break;
            }
        }
        if (!integrator_found){
            reb_simulation_warning(r,"Unknown integrator encountered in Simulationarchive. Reset function pointers manually.");
        }
    }

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
    r->ri_whfast512.recalculate_constants = 1;
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

