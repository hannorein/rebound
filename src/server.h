/**
 * @file    server.c
 * @brief   Opens a webserver to allow for platform independent visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. Part of the code is by Dave O'Hallaron, Carnegie Mellon (tiny.c).
 * 
 * @section LICENSE
 * Copyright (c) 2023 Hanno Rein, Dave O'Hallaron, Carnegie Mellon
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
#ifndef _SERVER_H
#define _SERVER_H
struct reb_server_data {
    struct reb_simulation* r;
    void* screenshot; // Screenshot data received by server (decoded)
    size_t N_screenshot; // Size of decoded screenshot data
    enum REB_STATUS status_before_screenshot;
    int port;
    int need_copy;
    int ready;
#ifdef SERVER
    int mutex_locked_by_integrate;  // Let's heartbeat find out if it is being called while the mutex is locked.
#ifdef _WIN32
    SOCKET socket;
    HANDLE mutex;          // Mutex to allow for copying
#else // _WIN32
    int socket;
    pthread_mutex_t mutex;          // Mutex to allow for copying
    pthread_t server_thread;
#endif // _WIN32
#endif // SERVER
};

#endif // _SERVER_H
