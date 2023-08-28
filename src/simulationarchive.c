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
    // reb_create_simulation sets simulationarchive_version to 3 by default.
    // This will break reading in old version.
    // Set to old version by default. Will be overwritten if new version was used.
    r->simulationarchive_version = 0;

    fseek(inf, 0, SEEK_SET);
    reb_input_fields(r, inf, warnings,NULL);

    // Done?
    if (snapshot==0) return;

    // Read SA snapshot
    if(fseek(inf, sa->offset64[snapshot], SEEK_SET)){
        *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
        reb_free_simulation(r);
        return;
    }
    if (r->simulationarchive_version<2){ 
        *warnings |= REB_INPUT_BINARY_ERROR_OLD;
        reb_free_simulation(r);
        return;
    }else{
        // Version 2 or higher
        reb_input_fields(r, inf, warnings,NULL);
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

void reb_read_simulationarchive_from_stream_with_messages(struct reb_simulationarchive* sa, struct reb_simulationarchive* sa_index, enum reb_input_binary_messages* warnings){
    // Assumes sa->inf is set to an open stream
    const int debug = 0;
    if (sa->inf==NULL){
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    
    // Get version
    fseek(sa->inf, 0, SEEK_SET);  
    struct reb_binary_field field = {0};
    sa->version = 0;
    double t0 = 0;
    
    // Cache descriptors
    struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
    struct reb_binary_field_descriptor fd_t = reb_binary_field_descriptor_for_name("t");
    struct reb_binary_field_descriptor fd_sa_version = reb_binary_field_descriptor_for_name("simulationarchive_version");
    struct reb_binary_field_descriptor fd_sa_size_snapshot = reb_binary_field_descriptor_for_name("simulationarchive_size_snapshot");
    struct reb_binary_field_descriptor fd_sa_size_first = reb_binary_field_descriptor_for_name("simulationarchive_size_first");
    struct reb_binary_field_descriptor fd_sa_auto_walltime = reb_binary_field_descriptor_for_name("simulationarchive_auto_walltime");
    struct reb_binary_field_descriptor fd_sa_auto_interval = reb_binary_field_descriptor_for_name("simulationarchive_auto_interval");
    struct reb_binary_field_descriptor fd_sa_auto_step = reb_binary_field_descriptor_for_name("simulationarchive_auto_step");
    struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");


    do{
        fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
        if (field.type == fd_header.type){
            long objects = 0;
            // Input header.
            const long bufsize = 64 - sizeof(struct reb_binary_field);
            char readbuf[bufsize], curvbuf[bufsize];
            const char* header = "REBOUND Binary File. Version: ";
            sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);

            objects += fread(readbuf,sizeof(char),bufsize,sa->inf);
            if (objects < 1){
                *warnings |= REB_INPUT_BINARY_WARNING_CORRUPTFILE;
            }else{
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                }
            }

        }else if (field.type == fd_t.type){
            fread(&t0, sizeof(double),1,sa->inf);
        }else if (field.type == fd_sa_version.type){
            fread(&(sa->version), sizeof(int),1,sa->inf);
        }else if (field.type == fd_sa_size_snapshot.type){
            fread(&(sa->size_snapshot), sizeof(long),1,sa->inf);
        }else if (field.type == fd_sa_size_first.type){
            fread(&(sa->size_first), sizeof(long),1,sa->inf);
        }else if (field.type == fd_sa_auto_walltime.type){
            fread(&(sa->auto_walltime), sizeof(double),1,sa->inf);
        }else if (field.type == fd_sa_auto_interval.type){
            fread(&(sa->auto_interval), sizeof(double),1,sa->inf);
        }else if (field.type == fd_sa_auto_step.type){
            fread(&(sa->auto_step), sizeof(unsigned long long),1,sa->inf);
        }else{
            fseek(sa->inf,field.size,SEEK_CUR);
        }
    }while(field.type!=fd_end.type);

    // Make index
    if (sa->version<2){
        // Version 1 no longer supported
        free(sa->filename);
        fclose(sa->inf);
        *warnings |= REB_INPUT_BINARY_ERROR_OLD;
        return;
    }else{
        // New version
        if (debug) printf("=============\n");
        if (debug) printf("SA Version: 2\n");
        if (sa_index == NULL){ // Need to construct offset index from file.
            long nblobsmax = 1024;
            sa->t = malloc(sizeof(double)*nblobsmax);
            sa->offset64 = malloc(sizeof(uint64_t)*nblobsmax);
            fseek(sa->inf, 0, SEEK_SET);  
            sa->nblobs = 0;
            int read_error = 0;
            struct reb_binary_field_descriptor fd_header = reb_binary_field_descriptor_for_name("header");
            struct reb_binary_field_descriptor fd_t = reb_binary_field_descriptor_for_name("t");
            struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
            for(long i=0;i<nblobsmax;i++){
                struct reb_binary_field field = {0};
                sa->offset64[i] = ftell(sa->inf);
                int blob_finished = 0;
                do{
                    size_t r1 = fread(&field,sizeof(struct reb_binary_field),1,sa->inf);
                    if (r1==1){
                        if (field.type == fd_header.type){
                            if (debug) printf("SA Field. type=HEADER\n");
                            int s1 = fseek(sa->inf,64 - sizeof(struct reb_binary_field),SEEK_CUR);
                            if (s1){
                                read_error = 1;
                            }
                        }else if (field.type == fd_t.type){
                            size_t r2 = fread(&(sa->t[i]), sizeof(double),1,sa->inf);
                            if (debug) printf("SA Field. type=TIME      value=%.10f\n",sa->t[1]);
                            if (r2!=1){
                                read_error = 1;
                            }
                        }else if (field.type == fd_end.type){
                            if (debug) printf("SA Field. type=END   =====\n");
                            blob_finished = 1;
                        }else{
                            int s2 = fseek(sa->inf,field.size,SEEK_CUR);
                            if (debug) printf("SA Field. type=%-6d    size=%llu\n",field.type,(unsigned long long)field.size);
                            if (s2){
                                read_error = 1;
                            }
                        }
                    }else{
                        read_error = 1;
                    }
                }while(blob_finished==0 && read_error==0);
                if (read_error){
                    if (debug) printf("SA Error. Error while reading current blob.\n");
                    // Error during reading. Current snapshot is corrupt.
                    break;
                }
                // Everything looks normal so far. Attempt to read next blob
                if (sa->version<3) { // will be removed in a future release
                    struct reb_simulationarchive_blob16 blob = {0};
                    size_t r3 = fread(&blob, sizeof(struct reb_simulationarchive_blob16), 1, sa->inf);
                    if (r3!=1){ // Next snapshot is definitly corrupted. Assume current might also be.
                        if (debug) printf("SA Error. Error while reading next blob.\n");
                        read_error = 1;
                        break;
                    }
                    if (i>0){
                        // Checking the offsets. Acts like a checksum.
                        if (((long)blob.offset_prev) + ((long)sizeof(struct reb_simulationarchive_blob16)) != ftell(sa->inf) - ((long)sa->offset64[i]) ){
                            // Offsets don't work. Next snapshot is definitly corrupted. Assume current one as well.
                            if (debug) printf("SA Error. Offset mismatch: %lu != %ld.\n",blob.offset_prev + sizeof(struct reb_simulationarchive_blob16), (long unsigned int)(ftell(sa->inf) - sa->offset64[i]) );
                            read_error = 1;
                            break;
                        }
                    }
                    // All tests passed. Accept current snapshot. Increase blob count.
                    sa->nblobs = i+1;
                    if (blob.offset_next==0){
                        // Last blob. 
                        if (debug) printf("SA Reached final blob.\n");
                        break;
                    }
                    if (i==nblobsmax-1){ // Increase 
                        nblobsmax += 1024;
                        sa->t = realloc(sa->t,sizeof(double)*nblobsmax);
                        sa->offset64 = realloc(sa->offset64,sizeof(uint64_t)*nblobsmax);
                    }
                }else{
                    struct reb_simulationarchive_blob blob = {0};
                    size_t r3 = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa->inf);
                    if (r3!=1){ // Next snapshot is definitly corrupted. Assume current might also be.
                        if (debug) printf("SA Error. Error while reading next blob.\n");
                        read_error = 1;
                        break;
                    }
                    if (i>0){
                        // Checking the offsets. Acts like a checksum.
                        if (((long)blob.offset_prev )+ ((long)sizeof(struct reb_simulationarchive_blob)) != ftell(sa->inf) - ((long)sa->offset64[i]) ){
                            // Offsets don't work. Next snapshot is definitly corrupted. Assume current one as well.
                            if (debug) printf("SA Error. Offset mismatch: %lu != %ld.\n",blob.offset_prev + sizeof(struct reb_simulationarchive_blob), (long unsigned int)(ftell(sa->inf) - sa->offset64[i]) );
                            read_error = 1;
                            break;
                        }
                    }
                    // All tests passed. Accept current snapshot. Increase blob count.
                    sa->nblobs = i+1;
                    if (blob.offset_next==0){
                        // Last blob. 
                        if (debug) printf("SA Reached final blob.\n");
                        break;
                    }
                    if (i==nblobsmax-1){ // Increase 
                        nblobsmax += 1024;
                        sa->t = realloc(sa->t,sizeof(double)*nblobsmax);
                        sa->offset64 = realloc(sa->offset64,sizeof(uint64_t)*nblobsmax);
                    }
                }
            }
            if (read_error){
                if (sa->nblobs>0){
                    *warnings |= REB_INPUT_BINARY_WARNING_CORRUPTFILE;
                }else{
                    fclose(sa->inf);
                    free(sa->filename);
                    free(sa->t);
                    free(sa->offset64);
                    free(sa);
                    *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
                    return;
                }
            }

        }else{ // reuse index from other SA
            // This is an optimzation for loading many large SAs.
            // It assumes the structure of this SA is *exactly* the same as in sa_index.
            // Unexpected behaviour if the shape is not the same.
            sa->nblobs = sa_index->nblobs;
            sa->t = malloc(sizeof(double)*sa->nblobs);
            sa->offset64 = malloc(sizeof(uint64_t)*sa->nblobs);
            fseek(sa->inf, 0, SEEK_SET);
            // No need to read the large file, just copying the index.
            memcpy(sa->offset64, sa_index->offset64, sizeof(uint64_t)*sa->nblobs);
            memcpy(sa->t, sa_index->t, sizeof(double)*sa->nblobs);
        }
    }
}

void reb_read_simulationarchive_with_messages(struct reb_simulationarchive* sa, const char* filename,  struct reb_simulationarchive* sa_index, enum reb_input_binary_messages* warnings){
    // Somewhat complicated calls for backwards compatability.
    sa->inf = fopen(filename,"r");
    sa->filename = malloc(strlen(filename)+1);
    strcpy(sa->filename,filename);
    reb_read_simulationarchive_from_stream_with_messages(sa, sa_index, warnings);
}

void reb_read_simulationarchive_from_buffer_with_messages(struct reb_simulationarchive* sa, char* buf, size_t size, struct reb_simulationarchive* sa_index, enum reb_input_binary_messages* warnings){
    // Somewhat complicated calls for backwards compatability.
    sa->inf = fmemopen(buf,size,"rb");
    sa->filename = NULL;
    reb_read_simulationarchive_from_stream_with_messages(sa, sa_index, warnings);
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
    free(sa->offset64);
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

void reb_simulationarchive_snapshot(struct reb_simulation* const r, const char* filename){
    if (r->simulationarchive_version<2){
        reb_error(r, "Writing SimulationArchives with version < 2 are no longer supported.\n");
        return;
    }
    if (filename==NULL) filename = r->simulationarchive_filename;
    struct stat buffer;
    if (stat(filename, &buffer) < 0){
        // File does not exist. Output binary.
        reb_output_binary(r,filename);
    }else{
        // File exists, append snapshot.
        struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
        if (r->simulationarchive_version<3){ // duplicate for working with old files. Will be removed in future release
                                             // Create buffer containing original binary file
            FILE* of = fopen(filename,"r+b");
            fseek(of, 64, SEEK_SET); // Header
            struct reb_binary_field field = {0};
            struct reb_simulationarchive_blob16 blob = {0};
            int bytesread;
            do{
                bytesread = fread(&field,sizeof(struct reb_binary_field),1,of);
                fseek(of, field.size, SEEK_CUR);
            }while(field.type!=fd_end.type && bytesread);
            long size_old = ftell(of);
            if (bytesread!=1){
                reb_warning(r, "SimulationArchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
                return;
            }

            bytesread = fread(&blob,sizeof(struct reb_simulationarchive_blob16),1,of);
            if (bytesread!=1){
                reb_warning(r, "SimulationArchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
                return;
            }
            int archive_contains_more_than_one_blob = 0;
            if (blob.offset_next>0){
                archive_contains_more_than_one_blob = 1;
            }

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

            int file_corrupt = 0;
            int seek_ok = fseek(of, -sizeof(struct reb_simulationarchive_blob16), SEEK_END);
            int blobs_read = fread(&blob, sizeof(struct reb_simulationarchive_blob16), 1, of);
            if (seek_ok !=0 || blobs_read != 1){ // cannot read blob
                file_corrupt = 1;
            }
            if ( (archive_contains_more_than_one_blob && blob.offset_prev <=0) || blob.offset_next != 0){ // blob contains unexpected data. Note: First blob is all zeros.
                file_corrupt = 1;
            }
            if (file_corrupt==0 && archive_contains_more_than_one_blob ){
                // Check if last two blobs are consistent.
                seek_ok = fseek(of, - sizeof(struct reb_simulationarchive_blob16) - sizeof(struct reb_binary_field), SEEK_CUR);  
                bytesread = fread(&field, sizeof(struct reb_binary_field), 1, of);
                if (seek_ok!=0 || bytesread!=1){
                    file_corrupt = 1;
                }
                if (field.type != fd_end.type || field.size !=0){
                    // expected an END field
                    file_corrupt = 1;
                }
                seek_ok = fseek(of, -blob.offset_prev - sizeof(struct reb_simulationarchive_blob16), SEEK_CUR);  
                struct reb_simulationarchive_blob16 blob2 = {0};
                blobs_read = fread(&blob2, sizeof(struct reb_simulationarchive_blob16), 1, of);
                if (seek_ok!=0 || blobs_read!=1 || blob2.offset_next != blob.offset_prev){
                    file_corrupt = 1;
                }
            }

            if (file_corrupt){
                // Find last valid snapshot to allow for restarting and appending to archives where last snapshot was cut off
                reb_warning(r, "SimulationArchive appears to be corrupted. REBOUND will attempt to fix it before appending more snapshots.\n");
                int seek_ok;
                seek_ok = fseek(of, size_old, SEEK_SET);
                long last_blob = size_old + sizeof(struct reb_simulationarchive_blob16);
                struct reb_binary_field_descriptor fd_end = reb_binary_field_descriptor_for_name("end");
                do
                {
                    seek_ok = fseek(of, -sizeof(struct reb_binary_field), SEEK_CUR);
                    if (seek_ok != 0){
                        break;
                    }
                    bytesread = fread(&field, sizeof(struct reb_binary_field), 1, of);
                    if (bytesread != 1 || field.type != fd_end.type){ // could be EOF or corrupt snapshot
                        break;
                    }
                    bytesread = fread(&blob, sizeof(struct reb_simulationarchive_blob16), 1, of);
                    if (bytesread != 1){
                        break;
                    }
                    last_blob = ftell(of);
                    if (blob.offset_next>0){
                        seek_ok = fseek(of, blob.offset_next, SEEK_CUR);
                    }else{
                        break;
                    }
                    if (seek_ok != 0){
                        break;
                    }
                } while(1);

                // To append diff, seek to last valid location (=EOF if all snapshots valid)
                fseek(of, last_blob, SEEK_SET);
            }else{
                // File is not corrupt. Start at end to save time.
                fseek(of, 0, SEEK_END);  
            }

            // Update blob info and Write diff to binary file
            fseek(of, -sizeof(struct reb_simulationarchive_blob16), SEEK_CUR);  
            fread(&blob, sizeof(struct reb_simulationarchive_blob16), 1, of);
            blob.offset_next = size_diff+sizeof(struct reb_binary_field);
            fseek(of, -sizeof(struct reb_simulationarchive_blob16), SEEK_CUR);  
            fwrite(&blob, sizeof(struct reb_simulationarchive_blob16), 1, of);
            fwrite(buf_diff, size_diff, 1, of); 
            field.type = fd_end.type;
            field.size = 0;
            fwrite(&field,sizeof(struct reb_binary_field), 1, of);
            blob.index++;
            blob.offset_prev = blob.offset_next;
            blob.offset_next = 0;
            fwrite(&blob, sizeof(struct reb_simulationarchive_blob16), 1, of);

            fclose(of);
            free(buf_new);
            free(buf_old);
            free(buf_diff);
        }else{ // Duplicate (version 3 of SimulationArchive. This is the part that will remain. Above duplicate will be removed in future release.
               // Create buffer containing original binary file
            FILE* of = fopen(filename,"r+b");
            fseek(of, 64, SEEK_SET); // Header
            struct reb_binary_field field = {0};
            struct reb_simulationarchive_blob blob = {0};
            int bytesread;
            do{
                bytesread = fread(&field,sizeof(struct reb_binary_field),1,of);
                fseek(of, field.size, SEEK_CUR);
            }while(field.type!=fd_end.type && bytesread);
            long size_old = ftell(of);
            if (bytesread!=1){
                reb_warning(r, "SimulationArchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
                return;
            }

            bytesread = fread(&blob,sizeof(struct reb_simulationarchive_blob),1,of);
            if (bytesread!=1){
                reb_warning(r, "SimulationArchive appears to be corrupted. A recovery attempt has failed. No snapshot has been saved.\n");
                return;
            }
            int archive_contains_more_than_one_blob = 0;
            if (blob.offset_next>0){
                archive_contains_more_than_one_blob = 1;
            }


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

            int file_corrupt = 0;
            int seek_ok = fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);
            int blobs_read = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
            if (seek_ok !=0 || blobs_read != 1){ // cannot read blob
                file_corrupt = 1;
            }
            if ( (archive_contains_more_than_one_blob && blob.offset_prev <=0) || blob.offset_next != 0){ // blob contains unexpected data. Note: First blob is all zeros.
                file_corrupt = 1;
            }
            if (file_corrupt==0 && archive_contains_more_than_one_blob ){
                // Check if last two blobs are consistent.
                seek_ok = fseek(of, - sizeof(struct reb_simulationarchive_blob) - sizeof(struct reb_binary_field), SEEK_CUR);  
                bytesread = fread(&field, sizeof(struct reb_binary_field), 1, of);
                if (seek_ok!=0 || bytesread!=1){
                    file_corrupt = 1;
                }
                if (field.type != fd_end.type || field.size !=0){
                    // expected an END field
                    file_corrupt = 1;
                }
                seek_ok = fseek(of, -blob.offset_prev - sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
                struct reb_simulationarchive_blob blob2 = {0};
                blobs_read = fread(&blob2, sizeof(struct reb_simulationarchive_blob), 1, of);
                if (seek_ok!=0 || blobs_read!=1 || blob2.offset_next != blob.offset_prev){
                    file_corrupt = 1;
                }
            }

            if (file_corrupt){
                // Find last valid snapshot to allow for restarting and appending to archives where last snapshot was cut off
                reb_warning(r, "SimulationArchive appears to be corrupted. REBOUND will attempt to fix it before appending more snapshots.\n");
                int seek_ok;
                seek_ok = fseek(of, size_old, SEEK_SET);
                long last_blob = size_old + sizeof(struct reb_simulationarchive_blob);
                do
                {
                    seek_ok = fseek(of, -sizeof(struct reb_binary_field), SEEK_CUR);
                    if (seek_ok != 0){
                        break;
                    }
                    bytesread = fread(&field, sizeof(struct reb_binary_field), 1, of);
                    if (bytesread != 1 || field.type != fd_end.type){ // could be EOF or corrupt snapshot
                        break;
                    }
                    bytesread = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
                    if (bytesread != 1){
                        break;
                    }
                    last_blob = ftell(of);
                    if (blob.offset_next>0){
                        seek_ok = fseek(of, blob.offset_next, SEEK_CUR);
                    }else{
                        break;
                    }
                    if (seek_ok != 0){
                        break;
                    }
                } while(1);

                // To append diff, seek to last valid location (=EOF if all snapshots valid)
                fseek(of, last_blob, SEEK_SET);
            }else{
                // File is not corrupt. Start at end to save time.
                fseek(of, 0, SEEK_END);  
            }

            // Update blob info and Write diff to binary file
            fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
            fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
            blob.offset_next = size_diff+sizeof(struct reb_binary_field);
            fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_CUR);  
            fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
            fwrite(buf_diff, size_diff, 1, of); 
            field.type = fd_end.type;
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
