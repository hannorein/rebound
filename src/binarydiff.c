/**
 * @file    binarydiff.c
 * @brief   Binary diff allows to compare binary snapshots.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
 * Copyright (c) 2018 Hanno Rein
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
#include "particle.h"
#include "rebound.h"
#include "tools.h"
#include "output.h"
#include "binarydiff.h"

// Wrapper for backwards compatibility
void reb_binary_diff(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep){
    // Ignores return value
    reb_binary_diff_with_options(buf1, size1, buf2, size2, bufp, sizep, 0);
}

int reb_binary_diff_with_options(char* buf1, size_t size1, char* buf2, size_t size2, char** bufp, size_t* sizep, int output_option){
    if (!buf1 || !buf2 || size1<64 || size2<64){
        printf("Cannot read input buffers.\n");
        return 0;
    }

    int are_different = 0;
    
    if (output_option==0){
        *bufp = NULL;
        *sizep = 0;
    }
    size_t allocatedsize = 0;

    // Header.
    if(memcmp(buf1,buf2,64)!=0){
        printf("Header in binary files are different.\n");
    }

    size_t pos1 = 64;
    size_t pos2 = 64;
    
    while(1){
        if (pos1+sizeof(struct reb_binary_field)>size1) break;
        struct reb_binary_field field1 = *(struct reb_binary_field*)(buf1+pos1);
        pos1 += sizeof(struct reb_binary_field);
        if (field1.type==REB_BINARY_FIELD_TYPE_END){
            break;
        }
        if (pos2+sizeof(struct reb_binary_field)>size2) pos2 = 64;
        struct reb_binary_field field2 = *(struct reb_binary_field*)(buf2+pos2);
        pos2 += sizeof(struct reb_binary_field);
        
        // Fields might not be in the same order.
        if (field1.type!=field2.type){
            // Will search for element in buf2, starting at beginning just past header
            // Note that we ignore all ADDITIONAL fields in buf2 that were not present in buf1 
            pos2 = 64;
            int notfound = 0; 
            while(1) {
                if (pos2+sizeof(struct reb_binary_field)>size2){
                    notfound = 1;
                    break;
                }
                field2 = *(struct reb_binary_field*)(buf2+pos2);
                pos2 += sizeof(struct reb_binary_field);
                if(field2.type==REB_BINARY_FIELD_TYPE_END){
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
                // Output field with size 0
                pos1 += field1.size; // For next search
                pos2 = 64; // For next search
                field1.size = 0;
                are_different = 1.;
                switch(output_option){
                    case 0:
                        reb_output_stream_write(bufp, &allocatedsize, sizep, &field1,sizeof(struct reb_binary_field));
                        break;
                    case 1:
                        printf("Field %d not in simulation 2.\n",field1.type);
                        break;
                    default:
                        break;
                }
                continue;
            }
        }
        // Can assume field1.type == field2.type from here on
        if (pos1+field1.size>size1) printf("Corrupt binary file buf1.\n");
        if (pos2+field2.size>size2) printf("Corrupt binary file buf2.\n");
        int fields_differ = 0;
        if (field1.size==field2.size){
            if (memcmp(buf1+pos1,buf2+pos2,field1.size)!=0){
                fields_differ = 1;
            }
        }else{
            fields_differ = 1;
        }
        if(fields_differ){
            if (field1.type!=REB_BINARY_FIELD_TYPE_WALLTIME){
                // Ignore the walltime field for the return value.
                // Typically we do not care about this field when comparing simulations.
                are_different = 1.;
            }
            switch(output_option){
                case 0:
                    reb_output_stream_write(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binary_field));
                    reb_output_stream_write(bufp, &allocatedsize, sizep, buf2+pos2,field2.size);
                    break;
                case 1:
                    printf("Field %d differs.\n",field1.type);
                    break;
                default:
                    break;
            }
        }
        pos1 += field1.size;
        pos2 += field2.size;
    }
    // Search for fields which are present in buf2 but not in buf1
    pos1 = 64;
    pos2 = 64;
    while(1){
        if (pos2+sizeof(struct reb_binary_field)>size2) break;
        struct reb_binary_field field2 = *(struct reb_binary_field*)(buf2+pos2);
        pos2 += sizeof(struct reb_binary_field);
        if (field2.type==REB_BINARY_FIELD_TYPE_END){
            break;
        }
        if (pos1+sizeof(struct reb_binary_field)>size1) pos1 = 64;
        struct reb_binary_field field1 = *(struct reb_binary_field*)(buf1+pos1);
        pos1 += sizeof(struct reb_binary_field);
        
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
            if (pos1+sizeof(struct reb_binary_field)>size1){
                notfound = 1;
                break;
            }
            field1 = *(struct reb_binary_field*)(buf1+pos1);
            pos1 += sizeof(struct reb_binary_field);
            if(field1.type==REB_BINARY_FIELD_TYPE_END){
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
        switch(output_option){
            case 0:
                reb_output_stream_write(bufp, &allocatedsize, sizep, &field2,sizeof(struct reb_binary_field));
                reb_output_stream_write(bufp, &allocatedsize, sizep, buf2+pos2,field2.size);
                break;
            case 1:
                printf("Field %d not in simulation 1.\n",field2.type);
                break;
            default:
                break;
        }
        pos1 = 64;
        pos2 += field2.size;
    }
    return are_different;
}
