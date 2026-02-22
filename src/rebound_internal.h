/**
 * @file    rebound_internal.h
 * @brief   Various functions used internally and implemented in rebound.c
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

#ifndef _REBOUND_INTERNAL_H
#define _REBOUND_INTERNAL_H

// Definitions for functions that we need to implement ourselves on Windows.
#ifdef _WIN32
typedef struct reb_timeval {
    int64_t tv_sec;
    int64_t tv_usec;
} reb_timeval;
int gettimeofday(struct reb_timeval * tp, struct timezone * tzp);
int asprintf(char **strp, const char *fmt, ...);
int rand_r (unsigned int *seed);
#include <io.h>
#else // Linux and MacOS
#define reb_timeval timeval
#include <sys/time.h>
#include <unistd.h>
#include <pthread.h>
#endif // _WIN32

// Githash should be provided as a command line argument to the compiler. 
// If not, use this dummy.
#ifndef GITHASH
#define GITHASH notavailable0000000000000000000000000001 
#endif

// Graceful global interrupt handler 
#include <signal.h>
extern volatile sig_atomic_t reb_sigint;  

// Global constants and variables
DLLEXPORT extern const char* reb_build_str;   ///< Date and time build string.
DLLEXPORT extern const char* reb_version_str; ///< Version string.
DLLEXPORT extern const char* reb_githash_str; ///< Current git hash.
DLLEXPORT extern const char* reb_logo[26];    ///< Logo of rebound. 
DLLEXPORT extern const unsigned char reb_favicon_png[]; /// < Favicon in PNG format.
DLLEXPORT extern const unsigned int reb_favicon_len;
DLLEXPORT extern const int reb_messages_max_length;
DLLEXPORT extern const int reb_messages_max_N;

// Only free memory in pointers of a simulation, but not the simulation itself.
DLLEXPORT void reb_simulation_free_pointers(struct reb_simulation* const r);

#endif // _REBOUND_INTERNAL_H
