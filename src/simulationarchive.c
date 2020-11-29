/**
 * @file 	simulationarchive.c
 * @brief 	Tools for creating and reading Simulation Archive binary files.
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdint.h>
#include "particle.h"
#include "rebound.h"
#include "binarydiff.h"
#include "output.h"
#include "tools.h"
#include "input.h"
#include "output.h"
#include "integrator_ias15.h"


void reb_create_simulation_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, long snapshot, enum reb_input_binary_messages* warnings){
    FILE* inf = sa->inf;
    if (inf == NULL){
        *warnings |= REB_INPUT_BINARY_ERROR_FILENOTOPEN;
        return;
    }
    if (snapshot<0) snapshot += sa->nblobs;
    if (snapshot>=sa->nblobs || snapshot<0){
        *warnings |= REB_INPUT_BINARY_ERROR_OUTOFRANGE;
        return;
    }
    
    // load original binary file
    reb_free_pointers(r);
    memset(r,0,sizeof(struct reb_simulation));
    reb_init_simulation(r);
    r->simulationarchive_filename = NULL;
    // reb_create_simulation sets simulationarchive_version to 2 by default.
    // This will break reading in old version.
    // Set to old version by default. Will be overwritten if new version was used.
    r->simulationarchive_version = 0;

    fseek(inf, 0, SEEK_SET);
    while(reb_input_field(r, inf, warnings,NULL)){ }

    // Done?
    if (snapshot==0) return;

    // Read SA snapshot
    if(fseek(inf, sa->offset[snapshot], SEEK_SET)){
        *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
        reb_free_simulation(r);
        return;
    }
    if (r->simulationarchive_version<2){ 
        fread(&(r->t),sizeof(double),1,inf);
        fread(&(r->walltime),sizeof(double),1,inf);
        if (r->simulationarchive_auto_interval){
            while (r->simulationarchive_next<=r->t){
                r->simulationarchive_next += r->simulationarchive_auto_interval;
            }
        }
        if (r->simulationarchive_auto_walltime){
            r->simulationarchive_next = r->walltime + r->simulationarchive_auto_interval;
        }
        switch (r->integrator){
            case REB_INTEGRATOR_JANUS:
                {
                    if (r->ri_janus.allocated_N<(unsigned int)r->N){
                        if (r->ri_janus.p_int){
                            free(r->ri_janus.p_int);
                        }
                        r->ri_janus.p_int= malloc(sizeof(struct reb_particle)*r->N);
                        r->ri_janus.allocated_N = r->N;
                    }
                    fread(r->ri_janus.p_int,sizeof(struct reb_particle_int)*r->N,1,inf);
                    reb_integrator_synchronize(r);  // get floating point coordinates 
                }
                break;
            case REB_INTEGRATOR_WHFAST:
            case REB_INTEGRATOR_SABA:
                {
                    // Recreate Jacobi arrrays
                    struct reb_particle* ps = r->particles;
                    if (r->ri_whfast.safe_mode==0){
                        // If same mode is off, store unsynchronized Jacobi coordinates
                        if (r->ri_whfast.allocated_N<(unsigned int)r->N){
                            if (r->ri_whfast.p_jh){
                                free(r->ri_whfast.p_jh);
                            }
                            r->ri_whfast.p_jh= malloc(sizeof(struct reb_particle)*r->N);
                            r->ri_whfast.allocated_N = r->N;
                        }
                        ps = r->ri_whfast.p_jh;
                    }
                    for(int i=0;i<r->N;i++){
                        fread(&(r->particles[i].m),sizeof(double),1,inf);
                        fread(&(ps[i].x),sizeof(double),1,inf);
                        fread(&(ps[i].y),sizeof(double),1,inf);
                        fread(&(ps[i].z),sizeof(double),1,inf);
                        fread(&(ps[i].vx),sizeof(double),1,inf);
                        fread(&(ps[i].vy),sizeof(double),1,inf);
                        fread(&(ps[i].vz),sizeof(double),1,inf);
                    }
                    if (r->ri_whfast.safe_mode==0){
                        // Assume we are not synchronized
                        r->ri_whfast.is_synchronized=0.;
                        // Recalculate total mass
                        double msum = r->particles[0].m;
                        for (int i=1;i<r->N;i++){
                            r->ri_whfast.p_jh[i].m = r->particles[i].m;
                            r->ri_whfast.p_jh[i].r = r->particles[i].r;
                            msum += r->particles[i].m;
                        }
                        r->ri_whfast.p_jh[0].m = msum;
                        r->ri_whfast.p_jh[0].r = r->particles[0].r;
                    }
                }
                break;
            case REB_INTEGRATOR_MERCURIUS:
                {
                    // Recreate heliocentric arrrays
                    struct reb_particle* ps = r->particles;
                    if (r->ri_mercurius.safe_mode==0){
                        // If same mode is off, store unsynchronized Jacobi coordinates
                        if (r->ri_whfast.allocated_N<(unsigned int)r->N){
                            if (r->ri_whfast.p_jh){
                                free(r->ri_whfast.p_jh);
                            }
                            r->ri_whfast.p_jh= malloc(sizeof(struct reb_particle)*r->N);
                            r->ri_whfast.allocated_N = r->N;
                        }
                        ps = r->ri_whfast.p_jh;
                    }
                    for(int i=0;i<r->N;i++){
                        fread(&(r->particles[i].m),sizeof(double),1,inf);
                        fread(&(ps[i].x),sizeof(double),1,inf);
                        fread(&(ps[i].y),sizeof(double),1,inf);
                        fread(&(ps[i].z),sizeof(double),1,inf);
                        fread(&(ps[i].vx),sizeof(double),1,inf);
                        fread(&(ps[i].vy),sizeof(double),1,inf);
                        fread(&(ps[i].vz),sizeof(double),1,inf);
                    }
                    if (r->ri_mercurius.dcrit){
                        free(r->ri_mercurius.dcrit);
                    }
                    r->ri_mercurius.dcrit = malloc(sizeof(double)*r->N);
                    r->ri_mercurius.dcrit_allocatedN = r->N;
                    fread(r->ri_mercurius.dcrit,sizeof(double),r->N,inf);
                    if (r->ri_mercurius.safe_mode==0){
                        // Assume we are not synchronized
                        r->ri_mercurius.is_synchronized=0.;
                        // Recalculate total mass
                        double msum = r->particles[0].m;
                        for (int i=1;i<r->N;i++){
                            r->ri_whfast.p_jh[i].m = r->particles[i].m;
                            r->ri_whfast.p_jh[i].r = r->particles[i].r;
                            msum += r->particles[i].m;
                        }
                        r->ri_whfast.p_jh[0].m = msum;
                        r->ri_whfast.p_jh[0].r = r->particles[0].r;
                    }
                }
                break;
            case REB_INTEGRATOR_IAS15:
                {
                    fread(&(r->dt),sizeof(double),1,inf);
                    fread(&(r->dt_last_done),sizeof(double),1,inf);
                    struct reb_particle* ps = r->particles;
                    for(int i=0;i<r->N;i++){
                        fread(&(ps[i].m),sizeof(double),1,inf);
                        fread(&(ps[i].x),sizeof(double),1,inf);
                        fread(&(ps[i].y),sizeof(double),1,inf);
                        fread(&(ps[i].z),sizeof(double),1,inf);
                        fread(&(ps[i].vx),sizeof(double),1,inf);
                        fread(&(ps[i].vy),sizeof(double),1,inf);
                        fread(&(ps[i].vz),sizeof(double),1,inf);
                    }
                    reb_integrator_ias15_alloc(r);
                    const int N3 = r->N*3;
                    reb_read_dp7(&(r->ri_ias15.b)  ,N3,inf,NULL);
                    reb_read_dp7(&(r->ri_ias15.csb),N3,inf,NULL);
                    reb_read_dp7(&(r->ri_ias15.e)  ,N3,inf,NULL);
                    reb_read_dp7(&(r->ri_ias15.br) ,N3,inf,NULL);
                    reb_read_dp7(&(r->ri_ias15.er) ,N3,inf,NULL);
                    fread((r->ri_ias15.csx),sizeof(double)*N3,1,inf);
                    fread((r->ri_ias15.csv),sizeof(double)*N3,1,inf);
                }
                break;
            default:
                *warnings |= REB_INPUT_BINARY_ERROR_INTEGRATOR;
                reb_free_simulation(r);
                return;
        }
    }else{
        // Version 2
        while(reb_input_field(r, inf, warnings,NULL)){ }
    }
    return;
}


struct reb_simulation* reb_create_simulation_from_simulationarchive(struct reb_simulationarchive* sa, long snapshot){
    if (sa==NULL) return NULL;
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_create_simulation();
    reb_create_simulation_from_simulationarchive_with_messages(r, sa, snapshot, &warnings);
    r = reb_input_process_warnings(r, warnings);
    return r; // might be null if error occured
}

void reb_read_simulationarchive_with_messages(struct reb_simulationarchive* sa, const char* filename,  struct reb_simulationarchive* sa_index, enum reb_input_binary_messages* warnings){
    sa->inf = fopen(filename,"r");
    if (sa->inf==NULL){
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    sa->filename = malloc(strlen(filename)+1);
    strcpy(sa->filename,filename);
    
    // Get version
    fseek(sa->inf, 0, SEEK_SET);  
    struct reb_binary_field field = {0};
    sa->version = 0;
    double t0 = 0;
    do{
        fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
        switch (field.type){
            case REB_BINARY_FIELD_TYPE_HEADER:
                //fseek(sa->inf,64 - sizeof(struct reb_binary_field),SEEK_CUR);
                {
                    long objects = 0;
                    // Input header.
                    const long bufsize = 64 - sizeof(struct reb_binary_field);
                    char readbuf[bufsize], curvbuf[bufsize];
                    const char* header = "REBOUND Binary File. Version: ";
                    sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);
                    
                    objects += fread(readbuf,sizeof(char),bufsize,sa->inf);
                    // Note: following compares version, but ignores githash.
                    if(strncmp(readbuf,curvbuf,bufsize)!=0){
                        *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                    }
                }
                break;
            case REB_BINARY_FIELD_TYPE_T:
                fread(&t0, sizeof(double),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SAVERSION:
                fread(&(sa->version), sizeof(int),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SASIZESNAPSHOT:
                fread(&(sa->size_snapshot), sizeof(long),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SASIZEFIRST:
                fread(&(sa->size_first), sizeof(long),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SAAUTOWALLTIME:
                fread(&(sa->auto_walltime), sizeof(double),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SAAUTOINTERVAL:
                fread(&(sa->auto_interval), sizeof(double),1,sa->inf);
                break;
            case REB_BINARY_FIELD_TYPE_SAAUTOSTEP:
                fread(&(sa->auto_step), sizeof(unsigned long long),1,sa->inf);
                break;
            default:
                fseek(sa->inf,field.size,SEEK_CUR);
                break;
        }
    }while(field.type!=REB_BINARY_FIELD_TYPE_END);

    // Make index
    if (sa->version<2){
        // Old version
        if (sa->size_first==-1 || sa->size_snapshot==-1){
            free(sa->filename);
            fclose(sa->inf);
            *warnings |= REB_INPUT_BINARY_ERROR_OUTOFRANGE;
            return;
        }
        fseek(sa->inf, 0, SEEK_END);  
        sa->nblobs = (ftell(sa->inf)-sa->size_first)/sa->size_snapshot+1; // +1 accounts for first binary 
        sa->t = malloc(sizeof(double)*sa->nblobs);
        sa->offset = malloc(sizeof(uint32_t)*sa->nblobs);
        sa->t[0] = t0;
        sa->offset[0] = 0;
        for(long i=1;i<sa->nblobs;i++){
            double offset = sa->size_first+(i-1)*sa->size_snapshot;
            fseek(sa->inf, offset, SEEK_SET);  
            fread(&(sa->t[i]),sizeof(double), 1, sa->inf);
            sa->offset[i] = offset;
        }
    }else{
        // New version
        struct reb_simulationarchive_blob blob = {0};
        fseek(sa->inf, -sizeof(struct reb_simulationarchive_blob), SEEK_END);  
        fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa->inf);
        sa->nblobs = blob.index+1;
        sa->t = malloc(sizeof(double)*sa->nblobs);
        sa->offset = malloc(sizeof(uint32_t)*sa->nblobs);
        fseek(sa->inf, 0, SEEK_SET);  
        if (sa_index){
            // This is an optimzation for loading many large SAs.
            // It assumes the structure of this SA is *exactly* the same as in sa_index.
            // Unexpected behaviour if the shape is not the same.
            if (sa->nblobs!=sa_index->nblobs) {
                fclose(sa->inf);
                free(sa->filename);
                free(sa->t);
                free(sa->offset);
                free(sa);
                *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
                return;
            }
            // No need to read the large file, just copying the index.
            memcpy(sa->offset, sa_index->offset, sizeof(uint32_t)*sa->nblobs);
            memcpy(sa->t, sa_index->t, sizeof(double)*sa->nblobs);
        }else{
            for(long i=0;i<sa->nblobs;i++){
                struct reb_binary_field field = {0};
                sa->offset[i] = ftell(sa->inf);
                do{
                    fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
                    switch (field.type){
                        case REB_BINARY_FIELD_TYPE_HEADER:
                            fseek(sa->inf,64 - sizeof(struct reb_binary_field),SEEK_CUR);
                            break;
                        case REB_BINARY_FIELD_TYPE_T:
                            fread(&(sa->t[i]), sizeof(double),1,sa->inf);
                            break;
                        default:
                            fseek(sa->inf,field.size,SEEK_CUR);
                            break;
                            
                    }
                }while(field.type!=REB_BINARY_FIELD_TYPE_END);
                fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa->inf);
                if (i!=blob.index) {
                    fclose(sa->inf);
                    free(sa->filename);
                    free(sa->t);
                    free(sa->offset);
                    free(sa);
                    *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
                    return;
                }
            }
        }
    }
}

struct reb_simulationarchive* reb_open_simulationarchive(const char* filename){
    struct reb_simulationarchive* sa = malloc(sizeof(struct reb_simulationarchive)); 
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    reb_read_simulationarchive_with_messages(sa, filename, NULL, &warnings);
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        // Don't output an error if file does not exist, just return NULL.
        free(sa);
        sa = NULL;
    }else{
        reb_input_process_warnings(NULL, warnings);
    }
    return sa;
}

void reb_close_simulationarchive(struct reb_simulationarchive* sa){
    reb_free_simulationarchive_pointers(sa);
    free(sa);
}

void reb_free_simulationarchive_pointers(struct reb_simulationarchive* sa){
    if (sa==NULL) return;
    if (sa->inf){
        fclose(sa->inf);
    }
    free(sa->filename);
    free(sa->t);
    free(sa->offset);
}
    
static int reb_simulationarchive_snapshotsize(struct reb_simulation* const r){
    int size_snapshot = 0;
    switch (r->integrator){
        case REB_INTEGRATOR_JANUS:
            size_snapshot = sizeof(double)*2+sizeof(struct reb_particle_int)*r->N;
            break;
        case REB_INTEGRATOR_WHFAST:
        case REB_INTEGRATOR_SABA:
            size_snapshot = sizeof(double)*2+sizeof(double)*7*r->N;
            break;
        case REB_INTEGRATOR_MERCURIUS:
            size_snapshot = sizeof(double)*2+sizeof(double)*8*r->N;
            break;
        case REB_INTEGRATOR_IAS15:
            size_snapshot =  sizeof(double)*4  // time, walltime, dt, dt_last_done
                             +sizeof(double)*3*r->N*5*7  // dp7 arrays
                             +sizeof(double)*7*r->N      // particle m, pos, vel
                             +sizeof(double)*3*r->N*2;   // csx, csv
            break;
        default:
            reb_error(r,"Simulation archive not implemented for this integrator.");
            break;
    }
    return size_snapshot;
}

void reb_simulationarchive_heartbeat(struct reb_simulation* const r){
    if (r->simulationarchive_filename!=NULL){
        int modes = 0;
        if (r->simulationarchive_auto_interval!=0) modes++;
        if (r->simulationarchive_auto_walltime!=0.) modes++;
        if (r->simulationarchive_auto_step!=0) modes++;
        if (modes>1){
            reb_error(r,"Only use one of simulationarchive_auto_interval, simulationarchive_auto_walltime, or simulationarchive_auto_step");
        }
        if (r->simulationarchive_auto_interval!=0.){
            const double sign = r->dt>0.?1.:-1;
            if (sign*r->simulationarchive_next <= sign*r->t){
                r->simulationarchive_next += sign*r->simulationarchive_auto_interval;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        }
        if (r->simulationarchive_auto_step!=0.){
            if (r->simulationarchive_next_step <= r->steps_done){
                r->simulationarchive_next_step += r->simulationarchive_auto_step;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        }
        if (r->simulationarchive_auto_walltime!=0.){
            if (r->simulationarchive_next <= r->walltime){
                r->simulationarchive_next += r->simulationarchive_auto_walltime;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        } 
    }
}
static inline void reb_save_dp7_old(struct reb_dp7* dp7, const int N3, FILE* of){
    fwrite(dp7->p0,sizeof(double),N3,of);
    fwrite(dp7->p1,sizeof(double),N3,of);
    fwrite(dp7->p2,sizeof(double),N3,of);
    fwrite(dp7->p3,sizeof(double),N3,of);
    fwrite(dp7->p4,sizeof(double),N3,of);
    fwrite(dp7->p5,sizeof(double),N3,of);
    fwrite(dp7->p6,sizeof(double),N3,of);
}

void reb_simulationarchive_snapshot(struct reb_simulation* const r, const char* filename){
    if (filename==NULL) filename = r->simulationarchive_filename;
    struct stat buffer;
    if (stat(filename, &buffer) < 0){
        // File does not exist. Output binary.
        if (r->simulationarchive_version<2){
            // Old version
            r->simulationarchive_size_snapshot = reb_simulationarchive_snapshotsize(r);
        }
        reb_output_binary(r,filename);
    }else{
        // File exists, append snapshot.
        if (r->simulationarchive_version<2){
            // Old version
            FILE* of = fopen(filename,"r+");
            fseek(of, 0, SEEK_END);
            fwrite(&(r->t),sizeof(double),1, of);
            fwrite(&(r->walltime),sizeof(double),1, of);
            switch (r->integrator){
                case REB_INTEGRATOR_JANUS:
                    {
                        fwrite(r->ri_janus.p_int,sizeof(struct reb_particle_int)*r->N,1,of);
                    }
                    break;
                case REB_INTEGRATOR_WHFAST:
                    {
                        struct reb_particle* ps = r->particles;
                        if (r->ri_whfast.safe_mode==0){
                            ps = r->ri_whfast.p_jh;
                        }
                        for(int i=0;i<r->N;i++){
                            fwrite(&(r->particles[i].m),sizeof(double),1,of);
                            fwrite(&(ps[i].x),sizeof(double),1,of);
                            fwrite(&(ps[i].y),sizeof(double),1,of);
                            fwrite(&(ps[i].z),sizeof(double),1,of);
                            fwrite(&(ps[i].vx),sizeof(double),1,of);
                            fwrite(&(ps[i].vy),sizeof(double),1,of);
                            fwrite(&(ps[i].vz),sizeof(double),1,of);
                        }
                    }
                    break;
                case REB_INTEGRATOR_MERCURIUS:
                    {
                        struct reb_particle* ps = r->particles;
                        if (r->ri_mercurius.safe_mode==0){
                            ps = r->ri_whfast.p_jh;
                        }
                        for(int i=0;i<r->N;i++){
                            fwrite(&(r->particles[i].m),sizeof(double),1,of);
                            fwrite(&(ps[i].x),sizeof(double),1,of);
                            fwrite(&(ps[i].y),sizeof(double),1,of);
                            fwrite(&(ps[i].z),sizeof(double),1,of);
                            fwrite(&(ps[i].vx),sizeof(double),1,of);
                            fwrite(&(ps[i].vy),sizeof(double),1,of);
                            fwrite(&(ps[i].vz),sizeof(double),1,of);
                        }
                        fwrite(r->ri_mercurius.dcrit,sizeof(double),r->N,of);
                    }
                    break;
                case REB_INTEGRATOR_IAS15:
                    {
                        fwrite(&(r->dt),sizeof(double),1,of);
                        fwrite(&(r->dt_last_done),sizeof(double),1,of);
                        struct reb_particle* ps = r->particles;
                        const int N3 = r->N*3;
                        for(int i=0;i<r->N;i++){
                            fwrite(&(ps[i].m),sizeof(double),1,of);
                            fwrite(&(ps[i].x),sizeof(double),1,of);
                            fwrite(&(ps[i].y),sizeof(double),1,of);
                            fwrite(&(ps[i].z),sizeof(double),1,of);
                            fwrite(&(ps[i].vx),sizeof(double),1,of);
                            fwrite(&(ps[i].vy),sizeof(double),1,of);
                            fwrite(&(ps[i].vz),sizeof(double),1,of);
                        }
                        reb_save_dp7_old(&(r->ri_ias15.b)  ,N3,of);
                        reb_save_dp7_old(&(r->ri_ias15.csb),N3,of);
                        reb_save_dp7_old(&(r->ri_ias15.e)  ,N3,of);
                        reb_save_dp7_old(&(r->ri_ias15.br) ,N3,of);
                        reb_save_dp7_old(&(r->ri_ias15.er) ,N3,of);
                        fwrite((r->ri_ias15.csx),sizeof(double)*N3,1,of);
                        fwrite((r->ri_ias15.csv),sizeof(double)*N3,1,of);
                    }
                    break;
                default:
                    reb_error(r,"Simulation archive not implemented for this integrator.");
                    break;
            }
            fclose(of);
        }else{
            // New version with incremental outputs

            // Create buffer containing original binary file
            FILE* of = fopen(filename,"r+b");
            fseek(of, 64, SEEK_SET); // Header
            struct reb_binary_field field;
            int bytesread;
            do{
                bytesread = fread(&field,sizeof(struct reb_binary_field),1,of);
                fseek(of, field.size, SEEK_CUR);
            }while(field.type!=REB_BINARY_FIELD_TYPE_END && bytesread);
            long size_old = ftell(of);
            char* buf_old = malloc(size_old);
            fseek(of, 0, SEEK_SET);  
            fread(buf_old, size_old,1,of);

            // Create buffer containing current binary file
            char* buf_new;
            size_t size_new;
            reb_output_binary_to_stream(r, &buf_new, &size_new);
            
            // Create buffer containing diff
            char* buf_diff;
            size_t size_diff;
            reb_binary_diff(buf_old, size_old, buf_new, size_new, &buf_diff, &size_diff);
            
            // Update blob info and Write diff to binary file
            struct reb_simulationarchive_blob blob = {0};
            fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);  
            fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
            blob.offset_next = size_diff+sizeof(struct reb_binary_field);
            fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);  
            fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
            fwrite(buf_diff, size_diff, 1, of); 
            field.type = REB_BINARY_FIELD_TYPE_END;
            field.size = 0;
            fwrite(&field,sizeof(struct reb_binary_field), 1, of);
            blob.index++;
            blob.offset_prev = blob.offset_next;
            blob.offset_next = 0;
            fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);

            fclose(of);
            free(buf_new);
            free(buf_old);
            free(buf_diff);
        }
    }
}

static int _reb_simulationarchive_automate_set_filename(struct reb_simulation* const r, const char* filename){
    if (r==NULL) return -1;
    if (filename==NULL){
        reb_error(r, "Filename missing.");
        return -1;
    }
    struct stat buffer;
    if (stat(filename, &buffer) == 0){
        reb_warning(r, "File in use for SimulationArchive already exists. Snapshots will be appended.");
    }
    free(r->simulationarchive_filename);
    r->simulationarchive_filename = malloc((strlen(filename)+1)*sizeof(char));
    strcpy(r->simulationarchive_filename, filename);
    return 0;
}

void reb_simulationarchive_automate_interval(struct reb_simulation* const r, const char* filename, double interval){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    if(r->simulationarchive_auto_interval != interval){
        // Only update simulationarchive_next if interval changed. 
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_interval = interval;
        r->simulationarchive_next = r->t;
    }
}

void reb_simulationarchive_automate_walltime(struct reb_simulation* const r, const char* filename, double walltime){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    // Note that this will create two snapshots if restarted.
    r->simulationarchive_auto_walltime = walltime;
    r->simulationarchive_next = r->walltime;
}

void reb_simulationarchive_automate_step(struct reb_simulation* const r, const char* filename, unsigned long long step){
    if(_reb_simulationarchive_automate_set_filename(r,filename)<0) return;
    if(r->simulationarchive_auto_step != step){
        // Only update simulationarchive_next if interval changed. 
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_step = step;
        r->simulationarchive_next_step = r->steps_done;
    }
}
