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
#include <stdint.h>
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "input.h"
#include "output.h"
#include "integrator_ias15.h"



int reb_simulationarchive_load_snapshot(struct reb_simulation* r, char* filename, long snapshot){
    if (access(filename, F_OK) == -1) return -1;
    if (!r) return -2;
    if (snapshot==0){
        // load original binary file
        enum reb_input_binary_messages warnings = 0;
        reb_create_simulation_from_binary_with_messages(r,filename,&warnings);
        if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
            reb_error(r,"Cannot read binary file. Check filename and file contents.");
        }
        return 0;
    }
    
    FILE* fd = fopen(filename,"r");
    int fseekret = 0;
    if (snapshot<0){
        // Find latest snapshot
        fseekret = fseek(fd,-r->simulationarchive_size_snapshot,SEEK_END);
    }else{
        fseekret = fseek(fd,r->simulationarchive_size_first + (snapshot-1)*r->simulationarchive_size_snapshot,SEEK_SET);
    }
    if (fseekret){
        // Seek didn't work.
        fclose(fd);
        return -3;
    }
    fread(&(r->t),sizeof(double),1,fd);
    fread(&(r->simulationarchive_walltime),sizeof(double),1,fd);
    gettimeofday(&r->simulationarchive_time,NULL);
    if (r->simulationarchive_interval){
        while (r->simulationarchive_next<=r->t){
            r->simulationarchive_next += r->simulationarchive_interval;
        }
    }
    switch (r->integrator){
        case REB_INTEGRATOR_WHFAST:
            {
                // Recreate Jacobi arrrays
                struct reb_particle* ps = r->particles;
                if (r->ri_whfast.safe_mode==0){
                    // If same mode is off, store unsynchronized Jacobi coordinates
                    if (r->ri_whfast.allocated_N<r->N){
                        if (r->ri_whfast.p_j){
                            free(r->ri_whfast.p_j);
                        }
                        if (r->ri_whfast.eta){
                            free(r->ri_whfast.eta);
                        }
                        r->ri_whfast.p_j= malloc(sizeof(struct reb_particle)*r->N);
                        r->ri_whfast.eta= malloc(sizeof(double)*r->N);
                        r->ri_whfast.allocated_N = r->N;
                    }
                    ps = r->ri_whfast.p_j;
                }
                for(int i=0;i<r->N;i++){
                    fread(&(r->particles[i].m),sizeof(double),1,fd);
                    fread(&(ps[i].x),sizeof(double),1,fd);
                    fread(&(ps[i].y),sizeof(double),1,fd);
                    fread(&(ps[i].z),sizeof(double),1,fd);
                    fread(&(ps[i].vx),sizeof(double),1,fd);
                    fread(&(ps[i].vy),sizeof(double),1,fd);
                    fread(&(ps[i].vz),sizeof(double),1,fd);
                }
                if (r->ri_whfast.safe_mode==0){
                    // Assume we are not synchronized
                    r->ri_whfast.is_synchronized=0.;
                    // Recalculate Jacobi masses
                    r->ri_whfast.eta[0] = r->particles[0].m;
                    r->ri_whfast.p_j[0].m = r->particles[0].m;
                    for (unsigned int i=1;i<r->N;i++){
                        r->ri_whfast.eta[i] = r->ri_whfast.eta[i-1] + r->particles[i].m;
                        r->ri_whfast.p_j[i].m = r->particles[i].m;
                    }
                }
            }
            break;
        case REB_INTEGRATOR_WHFASTHELIO:
            {
                // Recreate Heliocentric arrrays
                struct reb_particle* ps = r->particles;
                if (r->ri_whfasthelio.safe_mode==0){
                    // If same mode is off, store unsynchronized Heliocentric coordinates
                    if (r->ri_whfasthelio.allocated_N<r->N){
                        if (r->ri_whfasthelio.p_h){
                            free(r->ri_whfasthelio.p_h);
                        }
                        r->ri_whfasthelio.p_h= malloc(sizeof(struct reb_particle)*r->N);
                        r->ri_whfasthelio.allocated_N = r->N;
                    }
                    ps = r->ri_whfasthelio.p_h;
                }
                for(int i=0;i<r->N;i++){
                    fread(&(r->particles[i].m),sizeof(double),1,fd);
                    fread(&(ps[i].x),sizeof(double),1,fd);
                    fread(&(ps[i].y),sizeof(double),1,fd);
                    fread(&(ps[i].z),sizeof(double),1,fd);
                    fread(&(ps[i].vx),sizeof(double),1,fd);
                    fread(&(ps[i].vy),sizeof(double),1,fd);
                    fread(&(ps[i].vz),sizeof(double),1,fd);
                }
                if (r->ri_whfasthelio.safe_mode==0){
                    // Assume we are not synchronized
                    r->ri_whfasthelio.is_synchronized=0.;
                    // Recalculate Jacobi masses
                    r->ri_whfasthelio.p_h[0].m = r->particles[0].m;
                    for (unsigned int i=1;i<r->N;i++){
                        r->ri_whfasthelio.p_h[0].m += r->particles[i].m;
                        r->ri_whfasthelio.p_h[i].m = r->particles[i].m;
                    }
                }
            }
            break;
        case REB_INTEGRATOR_IAS15:
            {
                fread(&(r->dt),sizeof(double),1,fd);
                fread(&(r->dt_last_done),sizeof(double),1,fd);
                struct reb_particle* ps = r->particles;
                for(int i=0;i<r->N;i++){
                    fread(&(ps[i].m),sizeof(double),1,fd);
                    fread(&(ps[i].x),sizeof(double),1,fd);
                    fread(&(ps[i].y),sizeof(double),1,fd);
                    fread(&(ps[i].z),sizeof(double),1,fd);
                    fread(&(ps[i].vx),sizeof(double),1,fd);
                    fread(&(ps[i].vy),sizeof(double),1,fd);
                    fread(&(ps[i].vz),sizeof(double),1,fd);
                }
                reb_integrator_ias15_alloc(r);
                const int N3 = r->N*3;
                reb_read_dp7(&(r->ri_ias15.b)  ,N3,fd);
                reb_read_dp7(&(r->ri_ias15.csb),N3,fd);
                reb_read_dp7(&(r->ri_ias15.e)  ,N3,fd);
                reb_read_dp7(&(r->ri_ias15.br) ,N3,fd);
                reb_read_dp7(&(r->ri_ias15.er) ,N3,fd);
                fread((r->ri_ias15.csx),sizeof(double)*N3,1,fd);
                fread((r->ri_ias15.csv),sizeof(double)*N3,1,fd);
            }
            break;
        default:
            reb_error(r,"Simulation archive not implemented for this integrator.");
            break;
    }
    fclose(fd);

    return 0;
}

static int reb_simulationarchive_snapshotsize(struct reb_simulation* const r){
    int size_snapshot = 0;
    switch (r->integrator){
        case REB_INTEGRATOR_WHFAST:
        case REB_INTEGRATOR_WHFASTHELIO:
            size_snapshot = sizeof(double)*2+sizeof(double)*7*r->N;
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

long reb_simulationarchive_estimate_size(struct reb_simulation* const r, double tmax){
    if (r->simulationarchive_interval){
        long blobsize = reb_simulationarchive_snapshotsize(r);
        return blobsize*(long)ceil((tmax-r->t)/r->simulationarchive_interval);
    }else{
        reb_warning(r, "Variable simulationarchive_interval not set. Cannot estimate filesize.");
        return 0;
    }
}

struct reb_simulation* reb_create_simulation_from_simulationarchive(char* filename){
    if (access(filename, F_OK) == -1) return NULL;
    struct reb_simulation* r = reb_create_simulation_from_binary(filename);
    if (r){
        int ret = reb_simulationarchive_load_snapshot(r, filename, -1);
        if (ret){
            reb_warning(r,"Did not find any snapshots other than the initial one.");
        }
    }
    return r;
}

static void reb_simulationarchive_append(struct reb_simulation* r){
    FILE* of = fopen(r->simulationarchive_filename,"a");
    fwrite(&(r->t),sizeof(double),1, of);
    fwrite(&(r->simulationarchive_walltime),sizeof(double),1, of);
    switch (r->integrator){
        case REB_INTEGRATOR_WHFASTHELIO:
            {
                struct reb_particle* ps = r->particles;
                if (r->ri_whfasthelio.safe_mode==0){
                    ps = r->ri_whfasthelio.p_h;
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
        case REB_INTEGRATOR_WHFAST:
            {
                struct reb_particle* ps = r->particles;
                if (r->ri_whfast.safe_mode==0){
                    ps = r->ri_whfast.p_j;
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
                reb_save_dp7(&(r->ri_ias15.b)  ,N3,of);
                reb_save_dp7(&(r->ri_ias15.csb),N3,of);
                reb_save_dp7(&(r->ri_ias15.e)  ,N3,of);
                reb_save_dp7(&(r->ri_ias15.br) ,N3,of);
                reb_save_dp7(&(r->ri_ias15.er) ,N3,of);
                fwrite((r->ri_ias15.csx),sizeof(double)*N3,1,of);
                fwrite((r->ri_ias15.csv),sizeof(double)*N3,1,of);
            }
            break;
        default:
            reb_error(r,"Simulation archive not implemented for this integrator.");
            break;
    }
    fclose(of);
}

void reb_simulationarchive_heartbeat(struct reb_simulation* const r){
    if (r->simulationarchive_walltime==0){
        // First output
        r->simulationarchive_size_snapshot = reb_simulationarchive_snapshotsize(r);
        switch (r->gravity){
            case REB_GRAVITY_BASIC:
            case REB_GRAVITY_NONE:
                break;
            default:
                reb_error(r,"Simulation archive not implemented for this gravity module.");
                break;
        }
        r->simulationarchive_next = r->t + r->simulationarchive_interval;
        r->simulationarchive_walltime = 1e-300;
        gettimeofday(&r->simulationarchive_time,NULL);
        reb_output_binary(r,r->simulationarchive_filename);
    }else{
        // Appending outputs
        if (r->simulationarchive_interval){
            if (r->simulationarchive_next <= r->t){
                r->simulationarchive_next += r->simulationarchive_interval;
                
                struct timeval time_now;
                gettimeofday(&time_now,NULL);
                r->simulationarchive_walltime += time_now.tv_sec-r->simulationarchive_time.tv_sec+(time_now.tv_usec-r->simulationarchive_time.tv_usec)/1e6;
                r->simulationarchive_time = time_now;
                reb_simulationarchive_append(r);
            }
        }
        if (r->simulationarchive_interval_walltime){
            struct timeval time_now;
            gettimeofday(&time_now,NULL);
            double delta_walltime = time_now.tv_sec-r->simulationarchive_time.tv_sec+(time_now.tv_usec-r->simulationarchive_time.tv_usec)/1e6;
            if (delta_walltime >= r->simulationarchive_interval_walltime){
                r->simulationarchive_walltime += delta_walltime;
                r->simulationarchive_time = time_now;
                reb_simulationarchive_append(r);
            }
        }
    }
}
