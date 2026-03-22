/**
 * @file    display.h
 * @brief   Realtime OpenGL visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
 * Copyright (c) 2016 Hanno Rein, Shangfei Liu
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
#ifndef _DISPLAY_H
#define _DISPLAY_H

#ifdef OPENGL
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define GL_GLEXT_PROTOTYPES
#define EGL_EGLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#else // __EMSCRIPTEN__
#define GLFW_INCLUDE_NONE
#include "glad.h"
#endif // __EMSCRIPTEN__
#include <GLFW/glfw3.h>

struct reb_simulation;

struct reb_orbit_opengl {
    float x,y,z;
    float a, e, f;
    float omega, Omega, inc;
};
struct reb_vec4df {
    float x,y,z,r;
};

struct reb_display_data {
    struct reb_display_settings s;
    struct reb_simulation* r;
    struct reb_simulation* r_copy;
    void* screenshot; // Screenshot data to be sent to server
    struct reb_vec4df* particle_data;
    struct reb_orbit_opengl* orbit_data;
    uint64_t N_allocated;
    double mouse_x;
    double mouse_y;
    double retina;
    int take_one_screenshot;
#ifndef _WIN32
    int need_copy;
    pthread_mutex_t mutex;          // Mutex to allow for copying
    pthread_t compute_thread;
#endif // _WIN32
#ifdef __EMSCRIPTEN__
    int connection_status;
#endif
    uint64_t breadcrumb_last_steps_done;
    unsigned int breadcrumb_N_allocated;
    unsigned int breadcrumb_current_index;
    unsigned int mouse_action;      
    unsigned int key_mods;      
    unsigned int particle_buffer;
    unsigned int particle_buffer_current;
    unsigned int orbit_buffer;
    unsigned int orbit_buffer_current;
    void* window;
    struct {
        unsigned int texture;
        unsigned int program;
        unsigned int vao;
        unsigned int pos_location;
        unsigned int ypos_location;
        unsigned int scale_location;
        unsigned int aspect_location;
        unsigned int screen_aspect_location;
        unsigned int rotation_location;
        unsigned int texture_location;
        unsigned int charval_buffer;
    } shader_simplefont;
    struct {
        unsigned int program;
        unsigned int box_vao;
        unsigned int cross_vao;
        unsigned int ruler_vao;
        unsigned int mvp_location;
        unsigned int color_location;
    } shader_box;
    struct {
        unsigned int mvp_location;
        unsigned int color_location;
        unsigned int current_index_location;
        unsigned int breadcrumb_N_location;
        unsigned int N_location;
        unsigned int program;
        unsigned int particle_vao;
    } shader_point;
    struct {
        unsigned int mvp_location;
        unsigned int program;
        unsigned int particle_vao_current;
        unsigned int particle_vao;
    } shader_sphere;
    struct {
        unsigned int mvp_location;
        unsigned int current_index_location;
        unsigned int breadcrumb_N_location;
        unsigned int N_location;
        unsigned int vertex_count_location;
        unsigned int program;
        unsigned int particle_vao_current;
        unsigned int particle_vao;
        unsigned int vertex_count;
    } shader_orbit;
    struct {
        unsigned int mvp_location;
        unsigned int vertex_count_location;
        unsigned int program;
        unsigned int particle_vao_current;
        unsigned int vertex_count;
    } shader_plane;
};

void reb_display_init(struct reb_simulation* const r);

void reb_display_init_data(struct reb_simulation* const r);

#ifdef __EMSCRIPTEN__
void reb_display_keyboard(GLFWwindow* window, int key, int scancode, int action, int mods);
#endif // __EMSCRIPTEN__

#endif // OPENGL
#endif
