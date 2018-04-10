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
#include "binarydiff.h"


FILE* reb_binary_diff(FILE* f1, FILE* f2){
    
    if (!f1 || !f2){
        printf("Cannot read binary file.\n");
        return NULL;
    }

    long objects1 = 0;
    long objects2 = 0;
    // Input header.
    int bufN = 128;
    char* readbuf1 = malloc(sizeof(char)*bufN);
    char* readbuf2 = malloc(sizeof(char)*bufN);
    objects1 += fread(readbuf1,sizeof(char),64,f1);
    objects2 += fread(readbuf2,sizeof(char),64,f2);
    if(strcmp(readbuf1,readbuf2)!=0){
        printf("Header in binary files are different.\n");
    }
    
    int reading_fields1 = 1;
    int reading_fields2 = 1;
    while(reading_fields1==1 && reading_fields2==1){
        struct reb_binary_field field1;
        fread(&field1,sizeof(struct reb_binary_field),1,f1);
        if (field1.type==REB_BINARY_FIELD_TYPE_END){
            reading_fields1 = 0;
        }
        struct reb_binary_field field2;
        fread(&field2,sizeof(struct reb_binary_field),1,f2);
        if (field2.type==REB_BINARY_FIELD_TYPE_END){
            reading_fields2 = 0;
        }
        if (field1.type==field2.type){
            if (field1.size==field2.size){
                if (field1.size>bufN){
                    while (field1.size>bufN){
                        bufN *= 2;
                    }
                    readbuf1 = realloc(readbuf1, sizeof(char)*bufN);
                    readbuf2 = realloc(readbuf2, sizeof(char)*bufN);
                }
                fread(readbuf1, field1.size,1,f1);
                fread(readbuf2, field1.size,1,f2);
                
                if(memcmp(readbuf1,readbuf2,field1.size)!=0){
                    printf("Field %d differs.\n",field1.type);
                }
            }
        }
    }
    free(readbuf1);
    free(readbuf2);
    return NULL;
}
