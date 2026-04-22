/**
 * @file    rebound.c
 * @brief   Various REBOUND control structures and routine.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#define STDERR_FILENO 2
#else
#include <unistd.h>
#endif
#include "rebound.h"
#include "rebound_internal.h"
#ifdef OPENMP
#include <omp.h>
#endif
#define STRINGIFY(s) str(s)
#define str(s) #s
const size_t reb_messages_max_length = 1024;   // needs to be constant expression for array size
const size_t reb_messages_max_N = 10;
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "5.0.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.
const struct reb_integrator reb_integrator_none = {.name="none"};

// Allow python to discover all available integrators
const struct reb_integrator* reb_integrators_available[] = {
#define X(name) &reb_integrator_##name,
    REB_AVAILABLE_INTEGRATORS
#undef X
    NULL // Null terminated
};

void reb_exit(const char* const msg){
    // This function should also kill all children. 
    // Not implemented as pid is not easy to get to.
    // kill(pid, SIGKILL);
    reb_message(NULL, 0, REB_MESSAGE_TYPE_ERROR, msg);
    exit(EXIT_FAILURE);
}

void reb_message(char*** messages, int save_messages, enum REB_MESSAGE_TYPE type, const char* const msg){
    if (!save_messages || messages==NULL){
        fprintf(stderr,"\n");
        switch (type){
            case REB_MESSAGE_TYPE_INFO:
                fprintf(stderr, REB_STR_BOLD "REBOUND Message!" REB_STR_RESET);
                break;
            case REB_MESSAGE_TYPE_WARNING:
                if (isatty(STDERR_FILENO)) { // Is stderr a terminal?
                    fprintf(stderr, REB_STR_YELLOW_BOLD "Warning!" REB_STR_RESET);
                }else{
                    fprintf(stderr,"Warning!");
                }
                break;
            case REB_MESSAGE_TYPE_ERROR:
                if (isatty(STDERR_FILENO)) { // Is stderr a terminal?
                    fprintf(stderr, REB_STR_RED_BOLD "Error!" REB_STR_RESET);
                }else{
                    fprintf(stderr,"Error!");
                }
                break;
        }
        fprintf(stderr, " %s\n",msg);
    }else{
        // Note: not thread safe.
        if (*messages==NULL){
            *messages = calloc(reb_messages_max_N,sizeof(char*));
        }
        size_t n = 0;
        for (;n<reb_messages_max_N;n++){
            if ((*messages)[n]==NULL){
                break;
            }
        }
        if (n==reb_messages_max_N){
            free((*messages)[0]);
            for (size_t i=0;i<reb_messages_max_N-1;i++){
                (*messages)[i] = (*messages)[i+1];
            }
            (*messages)[reb_messages_max_N-1] = NULL;
            n= reb_messages_max_N-1;
        }
        (*messages)[n] = malloc(sizeof(char*)*reb_messages_max_length);
        // First character indicates type for python.
        (*messages)[n][0] = (char)type;
        strncpy((*messages)[n]+1, msg, reb_messages_max_length-2);
        (*messages)[n][reb_messages_max_length-1] = '\0';
    }
}

int reb_pop_message(char** messages, char* const buf){
    if (messages){
        char* w0 = messages[0];
        if (w0){
            for(size_t i=0;i<reb_messages_max_N-1;i++){
                messages[i] = messages[i+1];
            }
            messages[reb_messages_max_N-1] = NULL;
            strcpy(buf,w0);
            free(w0);
            return 1;
        }
    }
    return 0;
}

// Handles graceful shutdown. For example triggered by keyboard interrupts.
volatile sig_atomic_t reb_sigint;
void reb_sigint_handler(int signum) {
    if (signum == SIGINT){
        reb_sigint += 1;
    }
}

// Checks if floating point contractions are on. 
// If so, this will prevent unit tests from passing and bit-wise reproducibility will fail.
int reb_check_fp_contract(){
    double a = 1.2382309285234567;
    double b = 2.123478623874623234567;
    double c = 6.0284234234234567;

    double r1 = a*b+c;
    double ab = a*b;
    double r2 = ab+c;

    return r1 != r2;
}

// Wrapper to free pointers from python.
void reb_free(void* p){
    free(p);
}

#ifdef _WIN32

void PyInit_librebound() {};

// Source: https://codebrowser.dev/glibc/glibc/stdlib/rand_r.c.html 
int rand_r(unsigned int *seed) {
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next / 65536) % 2048;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    *seed = next;

    return result;
}


// Source: https://stackoverflow.com/a/40160038/115102
int vasprintf(char **strp, const char *fmt, va_list ap) {
    // _vscprintf tells you how big the buffer needs to be
    int len = _vscprintf(fmt, ap);
    if (len == -1) {
        return -1;
    }
    size_t size = (size_t)len + 1;
    char *str = malloc(size);
    if (!str) {
        return -1;
    }
    // _vsprintf_s is the "secure" version of vsprintf
    int r = vsprintf_s(str, len + 1, fmt, ap);
    if (r == -1) {
        free(str);
        return -1;
    }
    *strp = str;
    return r;
}
int asprintf(char **strp, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    int r = vasprintf(strp, fmt, ap);
    va_end(ap);
    return r;
}

// Source: https://stackoverflow.com/a/26085827/115102
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64
int gettimeofday(struct reb_timeval * tp, struct timezone * tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970
    static const uint64_t EPOCH = ((uint64_t) 116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime( &system_time );
    SystemTimeToFileTime( &system_time, &file_time );
    time =  ((uint64_t)file_time.dwLowDateTime )      ;
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec  = (int64_t) ((time - EPOCH) / 10000000L);
    tp->tv_usec = (int64_t) (system_time.wMilliseconds * 1000);
    return 0;
}

// Source: https://stackoverflow.com/a/17283549/115102
void usleep(__int64 usec){
    HANDLE timer;
    LARGE_INTEGER ft;

    ft.QuadPart = -(10*usec); // Convert to 100 nanosecond interval, negative value indicates relative time

    timer = CreateWaitableTimer(NULL, TRUE, NULL);
    SetWaitableTimer(timer, &ft, 0, NULL, NULL, 0);
    WaitForSingleObject(timer, INFINITE);
    CloseHandle(timer);
}
#endif // _WIN32

void reb_omp_set_num_threads(int num_threads){
#ifdef OPENMP
    omp_set_num_threads(num_threads);
#else // OPENMP
    (void)num_threads; // Unused.
    reb_message(NULL, 0, REB_MESSAGE_TYPE_ERROR, "Cannot set number of OpenMP threads. Recompile with OpenMP.");
#endif // OPENMP
}

const char* reb_logo[26] = {
    "          _                           _  ",
    "         | |                         | | ",
    " _ __ ___| |__   ___  _   _ _ __   __| | ",
    "| '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | ",
    "| | |  __/ |_) | (_) | |_| | | | | (_| | ",
    "|_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| ",
    "                                         ",
    "              `-:://::.`                 ",
    "          `/oshhoo+++oossso+:`           ",
    "       `/ssooys++++++ossssssyyo:`        ",
    "     `+do++oho+++osssso++++++++sy/`      ",
    "    :yoh+++ho++oys+++++++++++++++ss.     ",
    "   /y++hooyyooshooo+++++++++++++++oh-    ",
    "  -dsssdssdsssdssssssssssooo+++++++oh`   ",
    "  ho++ys+oy+++ho++++++++oosssssooo++so   ",
    " .d++oy++ys+++oh+++++++++++++++oosssod   ",
    " -h+oh+++yo++++oyo+++++++++++++++++oom   ",
    " `d+ho+++ys+++++oys++++++++++++++++++d   ",
    "  yys++++oy+++++++oys+++++++++++++++s+   ",
    "  .m++++++h+++++++++oys++++++++++++oy`   ",
    "   -yo++++ss++++++++++oyso++++++++oy.    ",
    "    .ss++++ho+++++++++++osys+++++yo`     ",
    "      :ss+++ho+++++++++++++osssss-       ",
    "        -ossoys++++++++++++osso.         ",
    "          `-/oyyyssosssyso+/.            ",
    "                ``....`                  ",
};

const unsigned char reb_favicon_png[] = {
    0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a, 0x00, 0x00, 0x00, 0x0d,
    0x49, 0x48, 0x44, 0x52, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x10,
    0x08, 0x06, 0x00, 0x00, 0x00, 0x1f, 0xf3, 0xff, 0x61, 0x00, 0x00, 0x00,
    0x01, 0x73, 0x52, 0x47, 0x42, 0x00, 0xae, 0xce, 0x1c, 0xe9, 0x00, 0x00,
    0x00, 0x44, 0x65, 0x58, 0x49, 0x66, 0x4d, 0x4d, 0x00, 0x2a, 0x00, 0x00,
    0x00, 0x08, 0x00, 0x01, 0x87, 0x69, 0x00, 0x04, 0x00, 0x00, 0x00, 0x01,
    0x00, 0x00, 0x00, 0x1a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xa0, 0x01,
    0x00, 0x03, 0x00, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0xa0, 0x02,
    0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x10, 0xa0, 0x03,
    0x00, 0x04, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00,
    0x00, 0x00, 0x34, 0x55, 0x71, 0xf2, 0x00, 0x00, 0x01, 0xaf, 0x49, 0x44,
    0x41, 0x54, 0x38, 0x11, 0x6d, 0xd3, 0x3f, 0x48, 0x55, 0x51, 0x1c, 0xc0,
    0xf1, 0xf7, 0xb4, 0x87, 0x29, 0x1a, 0xa4, 0x90, 0x6d, 0x89, 0xa2, 0x0e,
    0x49, 0x69, 0xbd, 0x10, 0x1a, 0x2a, 0x74, 0x14, 0x57, 0x77, 0x09, 0x5d,
    0x0c, 0x23, 0x52, 0x70, 0x11, 0x6c, 0x70, 0x73, 0x30, 0x78, 0x98, 0x36,
    0x94, 0x20, 0x51, 0xad, 0xd2, 0x5f, 0x88, 0x1c, 0x55, 0x5c, 0x24, 0x0c,
    0x27, 0x51, 0xb4, 0xa9, 0x10, 0x41, 0x6a, 0xe8, 0x9f, 0xf8, 0xfd, 0x5e,
    0xee, 0x79, 0x1e, 0xd4, 0x1f, 0x7c, 0xee, 0xf9, 0x77, 0xcf, 0xb9, 0xf7,
    0xfc, 0xee, 0xb9, 0xd9, 0xcc, 0xc9, 0xb8, 0x42, 0x57, 0x1f, 0x6e, 0xe1,
    0x07, 0xfe, 0xe2, 0x00, 0x1b, 0x78, 0x82, 0x75, 0x9c, 0x1a, 0x65, 0xf4,
    0x3e, 0x46, 0x01, 0x8d, 0xe8, 0xc5, 0x03, 0x18, 0x25, 0x18, 0xc6, 0x3e,
    0x26, 0xe1, 0xbd, 0x49, 0x38, 0x60, 0x94, 0xe3, 0x43, 0x5a, 0xce, 0x51,
    0x6e, 0x63, 0x16, 0xbe, 0xcd, 0x23, 0xbc, 0x83, 0x71, 0x09, 0xf3, 0x78,
    0x81, 0xb3, 0xc8, 0x64, 0xbd, 0x10, 0x53, 0x58, 0xc1, 0x6f, 0x34, 0xa1,
    0x15, 0x7b, 0xd8, 0x41, 0x1b, 0x7c, 0xd0, 0x2f, 0x78, 0xcf, 0x22, 0x9c,
    0xdc, 0x8d, 0xfb, 0x67, 0xb8, 0x5c, 0x83, 0x13, 0x9f, 0x23, 0x8e, 0x09,
    0x1a, 0x79, 0x98, 0x87, 0x87, 0xf8, 0x83, 0x1b, 0xe8, 0x80, 0xf9, 0x69,
    0xc1, 0x33, 0xdf, 0x60, 0x06, 0xde, 0xec, 0x42, 0x9d, 0xf8, 0x87, 0x2d,
    0x54, 0x63, 0x04, 0x97, 0xe1, 0xb8, 0x39, 0x58, 0x83, 0xe1, 0xbc, 0x21,
    0x34, 0xd8, 0xf8, 0x8a, 0xcf, 0x30, 0xf3, 0x35, 0xb8, 0x03, 0x33, 0x5f,
    0x8b, 0x10, 0xe7, 0xa9, 0xbc, 0x85, 0xdb, 0x8b, 0xe3, 0xbd, 0x8d, 0x65,
    0x54, 0x46, 0xbd, 0xd3, 0xd4, 0x07, 0x31, 0x16, 0xf5, 0x59, 0xbd, 0x00,
    0x13, 0x5d, 0x61, 0x23, 0x8d, 0x79, 0x93, 0xe3, 0x1e, 0x7f, 0xa6, 0x1d,
    0x17, 0x29, 0xcf, 0xa1, 0x80, 0xf6, 0xb4, 0x2f, 0x14, 0xdf, 0xa9, 0xb8,
    0x95, 0xd1, 0xd0, 0x41, 0x99, 0x75, 0x81, 0xd2, 0xa8, 0xa3, 0x87, 0xfa,
    0x4b, 0x78, 0x70, 0x5c, 0xb4, 0x0a, 0x71, 0x7c, 0xa2, 0xe1, 0xbe, 0x7d,
    0x90, 0x91, 0x73, 0x81, 0x6f, 0xa8, 0x83, 0x61, 0x76, 0x17, 0x92, 0x5a,
    0x26, 0xf3, 0x85, 0xf2, 0x6a, 0x5a, 0x8f, 0x0b, 0x93, 0x7e, 0x17, 0xf5,
    0xd8, 0x74, 0xc0, 0x6f, 0xee, 0xe9, 0x32, 0x4c, 0x54, 0x08, 0xbf, 0xf3,
    0xbd, 0xd0, 0x88, 0x4a, 0xdf, 0xd8, 0x5c, 0x38, 0xa7, 0xd5, 0x37, 0x58,
    0x45, 0x0e, 0xb7, 0x11, 0x6f, 0xc7, 0x43, 0x93, 0xc7, 0xf1, 0xf8, 0x4f,
    0xc7, 0x3e, 0xc2, 0xdc, 0x64, 0xdc, 0xb3, 0xfd, 0x0a, 0x1f, 0x93, 0xd6,
    0xd1, 0xc5, 0x27, 0xc5, 0x8b, 0x3a, 0xd2, 0x81, 0xd7, 0x48, 0x8e, 0x72,
    0x18, 0x74, 0xd5, 0x37, 0xb8, 0x8e, 0x2e, 0xf8, 0xe7, 0xed, 0xc2, 0x2f,
    0x72, 0x13, 0x4b, 0x68, 0xc6, 0x38, 0x1a, 0x30, 0x00, 0x4f, 0x6f, 0xf1,
    0x5f, 0xb0, 0x1e, 0xc2, 0x1f, 0xa8, 0x1f, 0x75, 0xf0, 0xc4, 0xb9, 0x35,
    0xcf, 0x8a, 0x07, 0xee, 0x29, 0xd6, 0x50, 0x8c, 0x43, 0x2d, 0xa4, 0x52,
    0x79, 0x8a, 0xe5, 0x13, 0x77, 0x00, 0x00, 0x00, 0x00, 0x49, 0x45, 0x4e,
    0x44, 0xae, 0x42, 0x60, 0x82
};
const unsigned int reb_favicon_len = 581;
