/**
 * @file    display.c
 * @brief   Realtime OpenGL visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. 
 * 
 * @section LICENSE
 * Copyright (c) 2015 Hanno Rein, Shangfei Liu
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
#include "rebound.h"
#include "rebound_internal.h"
#define DEG2RAD (M_PI/180.)
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#ifndef _WIN32
#include <pthread.h>
#endif // _WIN32
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include "display.h"
#include "binarydata.h"
#include "tools.h"
#include "particle.h"
#include "boundary.h"
#include "display.h"
#include "output.h"
#include "simulation.h"
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

static void reb_display_set_default_view(struct reb_simulation* const r, struct reb_display_settings* s){
    float  scale = 0.;
    // Need a scale for visualization
    if (r->root_size==-1){  
        scale = 0.;
        const struct reb_particle* p = r->particles;
        for (size_t i=0;i<r->N;i++){
            const double _r = sqrt(p[i].x*p[i].x+p[i].y*p[i].y+p[i].z*p[i].z);
            scale = MAX(scale, _r);
        }
        if(scale==0.){
            scale = 1.;
        }
        scale *= 1.1;
    }else{
        struct reb_vec3d boxsize = {
            .x = r->root_size*(double)r->N_root_x,
            .y = r->root_size*(double)r->N_root_y,
            .z = r->root_size*(double)r->N_root_z,
        };
        scale = MAX(boxsize.x, MAX(boxsize.y, boxsize.z));
    }

    struct reb_mat4df oldview = s->view;
    s->view = reb_mat4df_scale(reb_mat4df_identity(), 1./scale, 1./scale, 1./scale);
    if (oldview.m[1]==0. && oldview.m[2]==0. && oldview.m[4]==0. && oldview.m[6]==0.){
        struct reb_rotation rotation = {
            .ix = 1./sqrt(2.),
            .iy = 0.,
            .iz = 0.,
            .r = 1./sqrt(2.),
        };
        s->view = reb_mat4df_multiply(reb_rotation_to_mat4df(rotation), s->view);
    }else if (oldview.m[1]==0. && oldview.m[2]==0. && oldview.m[4]==0. && oldview.m[5]==0.){
        struct reb_rotation rotation = {
            .ix = 0.,
            .iy = -1./sqrt(2.),
            .iz = 0.,
            .r = 1./sqrt(2.),
        };
        s->view = reb_mat4df_multiply(reb_rotation_to_mat4df(rotation), s->view);
    }
}

void reb_display_settings_init(struct reb_simulation*r, struct reb_display_settings* s){
    size_t l1 = SIZE_MAX;
    size_t l2 = SIZE_MAX;
    reb_simulation_two_largest_particles(r, &l1, &l2);
    double largest_radius = 0;
    if (l1!=SIZE_MAX){
        largest_radius = r->particles[l1].r;
    }
    if (largest_radius > 0.0){
        s->spheres       = 1; 
    }else{
        s->spheres       = 0; 
    }
    s->pause             = 0; 
    s->multisample       = 1; 
    if (strcmp(r->integrator.name, "whfast")==0){
        s->wire          = 1; 
    }else{
        s->wire          = 0; 
    }
    s->breadcrumbs       = 0;
    s->onscreentext      = 1; 
    s->ghostboxes        = 0; 
    s->reference         = SIZE_MAX;
    s->view.m[1]=1; // this will make set_default_view show the xy plane
    reb_display_set_default_view(r, s);
}

void reb_simulation_add_orbit_colors(struct reb_simulation* r, struct reb_vec3df color) {
    struct reb_vec3df* new_color_list = calloc(r->N_orbit_colors + 1, sizeof(struct reb_vec3df));
    if(r->orbit_colors != NULL) {
        memcpy(new_color_list, r->orbit_colors, r->N_orbit_colors * sizeof(struct reb_vec3df));
        free(r->orbit_colors);
    }

    new_color_list[r->N_orbit_colors] = color;

    // Allocate and move the new_color_list to the orbit colors list
    r->orbit_colors = calloc(r->N_orbit_colors + 1, sizeof(struct reb_vec3df));
    memmove(r->orbit_colors, new_color_list, (r->N_orbit_colors + 1) * sizeof(struct reb_vec3df));
    r->N_orbit_colors++;
}

void reb_simulation_set_orbit_color(struct reb_simulation* r, struct reb_vec3df color, int i) {
    if(i < 0 || i >= r->N_orbit_colors) {
        printf("Invalid index for setting orbit color -> index %i is out of range of %zu colors.", i, r->N_orbit_colors);
        return;
    }

    r->orbit_colors[i] = color;
}

#ifdef OPENGL
#include "simplefont.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#include <emscripten/fetch.h>
// Need to use emscripten_ version of these functions because types are wrong otherwise
void emscripten_glVertexAttribDivisor(GLuint index, GLuint divisor);
void emscripten_glDrawArraysInstanced(GLenum mode, GLint first, GLsizei count, GLsizei instancecount);
void emscripten_glDrawElementsInstanced(GLenum mode, GLint first, GLsizei count, const void* indices, GLsizei instancecount);
#define reb_glVertexAttribDivisor emscripten_glVertexAttribDivisor
#define reb_glDrawArraysInstanced emscripten_glDrawArraysInstanced
#define reb_glDrawElementsInstanced emscripten_glDrawElementsInstanced

EM_JS(void, reb_overlay_update, (const char* text, int status), {
        var overlaytext = document.getElementById("overlaytext");
        if (overlaytext){
        overlaytext.innerHTML = UTF8ToString(text);
        }
        var overlay = document.getElementById("overlay");
        if (overlay){
        if (status==-3){ // Pause
        overlay.style.backgroundColor = "rgba(100.0, 100.0, 0.0, 0.5)";
        }else if (status==0 || status==5){ // Finished.
        overlay.style.backgroundColor = "rgba(0.0, 255.0, 0.0, 0.5)";
        }else if (status==10){ // Connection error.
        overlay.style.backgroundColor = "rgba(255.0, 0.0, 0.0, 0.5)";
        }else{
        overlay.style.backgroundColor = "rgba(0, 0, 0, 0.5)";
        }
        }
        });
EM_JS(int, reb_overlay_help_show, (int show), {
        var overlaytoggle = document.getElementById("overlay-toggle");
        if (overlaytoggle){
        if (overlaytoggle.innerHTML == "1"){
        overlaytoggle.innerHTML = "";
        show = !show;
        }
        }
        var overlayhelp = document.getElementById("overlay-help");
        if (show){
        overlayhelp.style.display = "block";
        }else{
        overlayhelp.style.display = "none";
        }
        return show;
        });

#else
#define reb_glVertexAttribDivisor glVertexAttribDivisor
#define reb_glDrawArraysInstanced glDrawArraysInstanced
#define reb_glDrawElementsInstanced glDrawElementsInstanced
#endif

void reb_render_frame(void* p);

static const char* onscreenhelp[] = { 
    "REBOUND mouse and keyboard commands",
    "----------------------------------------------------",
    " To rotate the view, simply drag the simulation",
    " with the mouse. To zoom in, press the shift key ",
    " and then drag the simulation with the mouse.",
    "----------------------------------------------------",
    " h       | Show/hide this page",
    " q       | Quit simulation",
    " (space) | Pause simulation",
    " (ar dwn)| Perform one single time step",
    " (pg dwn)| Perform 50 time steps",
#ifdef __EMSCRIPTEN__
    " e       | Take screenshot and export as png file",
#else // __EMSCRIPTEN__
    " e       | Take screenshot and export as tga file",
#endif // __EMSCRIPTEN__
    " d       | Pause real-time visualization", 
    "         | (the simulation continues)",
    " r       | Reset view. Press multiple times to",
    "         | change orientation",
    " x/X     | Move to a coordinate system centered",
    "         | on a particle (note: does not work if", 
    "         | particle array is resorted)",
    " t       | Show/hide logo, time, timestep, number",
    "         | of particles, and scale",
    " s       | Toggle points/spheres/points+spheres/none",
    " g       | Toggle ghost boxes",
    " m       | Toggle multisampling",
    " w       | Toggle orbit mode (none/wire/plane)",
    " i / o   | Increase / decrease number of breadcrumbs",
    " c       | Clear breadcrumb data",
    " b       | Cycle through thick wire views",
    " p       | Toggle glow rendering",
    "----------------------------------------------------"
};


static int convertLine(const char* in, float* out){
    int j = 0;
    while(in[j]!=0&&in[j]!=10){  // end on new line or \0
        out[j*2+0] = (float)((((int)in[j]))%16);
        out[j*2+1] = (float)((((int)in[j]))/16);
        j++;
    }
    return j;
}

static unsigned int compileShader(int* shader, int type, const char* source){
    GLint status;

    *shader = glCreateShader(type);
    glShaderSource(*shader, 1, &source, NULL);
    glCompileShader(*shader);

    GLint logLength;
    glGetShaderiv(*shader, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0) {
        GLchar *log = (GLchar *)malloc(logLength);
        glGetShaderInfoLog(*shader, logLength, &logLength, log);
        printf("\n\n%s\n\n",log);
        free(log);
    }

    glGetShaderiv(*shader, GL_COMPILE_STATUS, &status);
    if (status == 0) {
        glDeleteShader(*shader);
        return 0;
    }

    return 1;
}

static unsigned int linkProgram(int prog){
    GLint status;
    glLinkProgram(prog);

    GLint logLength;
    glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 1) { // 0 should work - bug in emscripten?
        GLchar *log = (GLchar *)malloc(logLength);
        glGetProgramInfoLog(prog, logLength, &logLength, log);
        printf("\n\n%s\n\n",log);
        free(log);
    }

    glGetProgramiv(prog, GL_LINK_STATUS, &status);
    if (status == 0) {
        return 0;
    }

    return 1;
}

static int loadShader(const char* vert_source, const char* frag_source){
    GLint vertShader, fragShader;

    GLuint _program = glCreateProgram();

    if(!compileShader(&vertShader, GL_VERTEX_SHADER, vert_source)) {
        printf("Failed to compile vertex shader.\n");
        return -1;
    }

    if (!compileShader(&fragShader, GL_FRAGMENT_SHADER, frag_source)){
        printf("Failed to compile fragment shader.\n");
        return -1;
    }

    glAttachShader(_program, vertShader);
    glAttachShader(_program, fragShader);

    if (!linkProgram(_program)) {
        printf("Failed to link shader.\n");
        return -1;
    }

    if (vertShader) {
        glDetachShader(_program, vertShader);
        glDeleteShader(vertShader);
    }
    if (fragShader) {
        glDetachShader(_program, fragShader);
        glDeleteShader(fragShader);
    }
    return _program;
}

static void reb_glfw_error_callback(int error, const char* description){
    fprintf(stderr, "GLFW Error (%d): %s\n", error, description);
}

static void reb_display_scroll(GLFWwindow* window, double xoffset, double yoffset){
    (void) xoffset;
    (void) yoffset;
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_scroll\n");
        return;
    }
    float scale = 1.-yoffset/100.;
    data->s.view = reb_mat4df_scale(data->s.view, scale, scale, scale);
}
static void reb_display_mouse_button(GLFWwindow* window, int button, int action, int mods){
    (void)button;
    (void)mods;
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_mouse_button\n");
        return;
    }
    data->mouse_action = action;
}

static void reb_display_resize(GLFWwindow* window, int x, int y){
    (void)x;
    (void)y;
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_resize\n");
        return;
    }
}

static void reb_display_cursor(GLFWwindow* window, double x, double y){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_cursor\n");
        return;
    }
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    if (data->mouse_action==GLFW_RELEASE){
        // No button pressed
        data->mouse_x = INFINITY;
        data->mouse_y = INFINITY;
        return;
    }
    if (isinf(data->mouse_x)){
        // New drag event
        data->mouse_x = x;
        data->mouse_y = y;
        return;
    }
    if (data->mouse_action==GLFW_PRESS){
        if ((data->key_mods&GLFW_MOD_SHIFT)==0){
            // Drag 
            float dx = 3.*(x-data->mouse_x)/width;
            float dy = 3.*(y-data->mouse_y)/height;
            struct reb_rotation rot_dy = {.ix=sin(dy), .r=cos(dy)};
            struct reb_rotation rot_dx = {.iy=sin(dx), .r=cos(dx)};
            data->s.view = reb_mat4df_multiply(reb_rotation_to_mat4df(rot_dy), data->s.view);
            data->s.view = reb_mat4df_multiply(reb_rotation_to_mat4df(rot_dx), data->s.view);
        }else{
            // Zoom
            float ix = data->mouse_x/width-0.5;
            float iy = data->mouse_y/height-0.5;
            float ir = sqrt(ix*ix + iy*iy);
            float nx = x/width-0.5;
            float ny = y/height-0.5;
            float nr = sqrt(nx*nx + ny*ny);
            data->s.view = reb_mat4df_scale(data->s.view, nr/ir, nr/ir, nr/ir);
        }
        data->mouse_x = x;
        data->mouse_y = y;
        return;
    }
}

#define xstr(s) ystr(s)
#define ystr(s) #s
static void reb_display_clear_particle_data(struct reb_display_data* data){
    int N = data->N_allocated;
    int N_hist = data->breadcrumb_N_allocated;
    if (data->particle_data){
        float n = NAN;
        for (int i=0; i<N; i++){
            data->particle_data[i].x = n;
            data->particle_data[i].y = n;
            data->particle_data[i].z = n;
            data->particle_data[i].r = n;
        }
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        for (int i=0; i<N_hist; i++){
            glBufferSubData(GL_ARRAY_BUFFER, i*N*sizeof(struct reb_vec4df), N*sizeof(struct reb_vec4df), data->particle_data);
        }
    }
    if (data->orbit_data){
        float n = NAN;
        for (int i=0; i<N; i++){
            data->orbit_data[i].x = n; // enough to not render
            data->orbit_data[i].y = n;
            data->orbit_data[i].z = n;
        }
        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
        for (int i=0; i<N_hist; i++){
            glBufferSubData(GL_ARRAY_BUFFER, i*(N-1)*sizeof(struct reb_orbit_opengl), (N-1)*sizeof(struct reb_orbit_opengl), data->orbit_data);
        }
    }
}

void reb_display_keyboard(GLFWwindow* window, int key, int scancode, int action, int mods){
    (void)scancode;
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_keyboard\n");
        return;
    }
    if (!data->r){
        printf("Error accessing data->r in reb_display_keyboard\n");
        return;
    }
    // User defined keys:
    int skip_default_keys = 0;
    if (data->r->key_callback){
        skip_default_keys = data->r->key_callback(data->r, key);
    } 
    if (skip_default_keys){
        return;
    }
    // Default keys:
    data->key_mods = mods;
    if (action==GLFW_PRESS){
        switch(key){
            case 'H':
                data->s.onscreenhelp = !data->s.onscreenhelp;
                break;
            case 'Q':
                data->r->status = REB_STATUS_USER;
                break;
            case ' ':
                if (data->r->status == REB_STATUS_PAUSED){
                    printf("Resume.\n");
                    data->r->status = REB_STATUS_RUNNING;
                }else if (data->r->status == REB_STATUS_RUNNING || data->r->status == REB_STATUS_LAST_STEP){
                    printf("Pause.\n");
                    data->r->status = REB_STATUS_PAUSED;
                }
                break;
            case 'S':
                data->s.spheres = (data->s.spheres+1)%4;
                break;
            case 'G':
                data->s.ghostboxes = !data->s.ghostboxes;
                break;
            case 'M':
                data->s.multisample = !data->s.multisample;
                if (data->s.multisample){
                    glEnable(GL_MULTISAMPLE); 
                }else{
                    glDisable(GL_MULTISAMPLE); 
                }
                break;
            case 'R':
                data->s.reference     = SIZE_MAX;
                reb_display_set_default_view(data->r, &data->s);
                break;
            case 'D':
                data->s.pause = !data->s.pause;
                break;
            case 'W':
                data->s.wire = (data->s.wire+1)%3;
                break;
            case 'C':
                reb_display_clear_particle_data(data);
                break;
            case 'E':
                data->take_one_screenshot = 1;
                break;
            case 'I':
                data->s.breadcrumbs = MAX(1,data->s.breadcrumbs*2);
                break;
            case 'O':
                data->s.breadcrumbs = MAX(0, data->s.breadcrumbs/2) ;
                data->breadcrumb_current_index = 0; // prevent bad memory access after rescale 
                break;
            case 'T':
                data->s.onscreentext = !data->s.onscreentext;
                break;
            case 'X': 
                if (mods!=GLFW_MOD_SHIFT){
                    data->s.reference++;
                    if (data->s.reference>=data->r->N) data->s.reference = data->r->N-1;
                }else{
                    data->s.reference--;
                    if (data->s.reference==SIZE_MAX-1) data->s.reference = SIZE_MAX;
                }
                if (data->s.reference == SIZE_MAX){
                    printf("Reference particle: NONE.\n");
                }else{
                    printf("Reference particle: %zu.\n",data->s.reference);
                }
                break;
            case 264: // arrow down
                if (data->r->status == REB_STATUS_PAUSED){
                    data->r->status = REB_STATUS_SINGLE_STEP;
                    printf("Step.\n");
                }
                break;
            case 267: // page down
                if (data->r->status == REB_STATUS_PAUSED){
                    data->r->status = REB_STATUS_SINGLE_STEP - 50;
                    printf("50 steps.\n");
                }
                break;
            
            case 'B':
                if(data->s.drawwiresthick) {
                    int last_focused_section = data->shader_thick_orbit.focused_orbit_section_to_use;
                    data->shader_thick_orbit.focused_orbit_section_to_use = (last_focused_section + 1) % 3;

                    if(last_focused_section == 2) {
                        data->s.drawwiresthick = 0;
                        data->shader_thick_orbit.focused_orbit_section_to_use = 0;
                    }
                }
                else {
                    data->s.drawwiresthick = 1;
                    data->shader_thick_orbit.focused_orbit_section_to_use = 0;
                }
                
                break;
            
            case 'P':
                data->s.drawglow = (data->s.drawglow + 1) % 2;
                if(!data->reb_post_processing.post_processing_enabled) {
                    printf("Post processing is disabled due to error.\n");
                    data->s.drawglow = 0;
                }
                break;
        }
    }
}

// Actual rendering
// Makes a copy of the simulation first.
void reb_render_frame(void* p){
    struct reb_display_data* data = (struct reb_display_data*)p;
    struct reb_simulation* r = data->r;
    if (!data){
        printf("reb_display_data undefined in reb_render_frame().\n");
        return;
    }
    int width, height;
#ifdef __EMSCRIPTEN__
    // Need to query canvas size using JS, set window size, then read framebuffer size.
    width = EM_ASM_INT({
            return document.getElementById("canvas").scrollWidth;
            });
    height = EM_ASM_INT({
            return document.getElementById("canvas").scrollHeight;
            });
#endif
    int cwidth, cheight;
    glfwGetWindowSize(data->window, &cwidth, &cheight);
#ifdef __EMSCRIPTEN__
    if (cwidth!=width || cheight !=height){
        glfwSetWindowSize(data->window, width, height);
    }
#endif
    glfwGetFramebufferSize(data->window, &width, &height);

    // Check if we have a retina display
    data->retina = (double)width/(double)cwidth;

    struct reb_simulation* r_copy = r->display_data->r_copy;
    if (!r_copy){
        data->r_copy = reb_simulation_create();
        r_copy = data->r_copy;
    }

    // lock mutex for update
    data->need_copy = 1;
    int wait_count = 0;
    const int wait_count_max = 10;
    int ret_try = EBUSY;
    while (wait_count<wait_count_max && ret_try){// wait at most one frame.
        ret_try = pthread_mutex_trylock(&data->mutex);
        if (ret_try){ // not locked
            usleep(1./120.*1e6/wait_count_max); 
            wait_count++;
        }
    }

    if (!ret_try){
        // Copy if lock obtained. Otherwise use old data.
        enum REB_BINARYDATA_ERROR_CODE warnings = REB_BINARYDATA_WARNING_NONE;
        reb_simulation_copy_with_messages(data->r_copy,r,&warnings);
        data->need_copy = 0;
        pthread_mutex_unlock(&(data->mutex));  
    }

    if (r_copy->display_settings){
        // User provided settings server-side. Will overwrite our own.
        data->s = *r_copy->display_settings;
    }

    // prepare data (incl orbit calculation)
    size_t N = r_copy->N;

    if (N > data->N_allocated || data->s.breadcrumbs+1 != data->breadcrumb_N_allocated){
        data->N_allocated = N;
        data->breadcrumb_N_allocated = data->s.breadcrumbs+1;

        data->particle_data = realloc(data->particle_data, data->N_allocated*sizeof(struct reb_vec4df));
        data->orbit_data = realloc(data->orbit_data, data->N_allocated*sizeof(struct reb_orbit_opengl));

        // Resize memory if needed
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        glBufferData(GL_ARRAY_BUFFER, data->breadcrumb_N_allocated*data->N_allocated*sizeof(struct reb_vec4df), NULL, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer_current);
        glBufferData(GL_ARRAY_BUFFER, data->N_allocated*sizeof(struct reb_vec4df), NULL, GL_STATIC_DRAW);

        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
        glBufferData(GL_ARRAY_BUFFER, data->breadcrumb_N_allocated*data->N_allocated*sizeof(struct reb_orbit_opengl), NULL, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer_current);
        glBufferData(GL_ARRAY_BUFFER, data->N_allocated*sizeof(struct reb_orbit_opengl), NULL, GL_STATIC_DRAW);

        reb_display_clear_particle_data(data);

        glBindBuffer(GL_ARRAY_BUFFER, data->shader_thick_orbit.thick_orbit_cbo);
        struct reb_vec3df colors_to_allocate[N];
        for(int i = 0; i < N; i++) {
            if(i < r_copy->N_orbit_colors) {
                struct reb_vec3df curr_color = r_copy->orbit_colors[i];
                colors_to_allocate[i].x = curr_color.x;
                colors_to_allocate[i].y = curr_color.y;
                colors_to_allocate[i].z = curr_color.z;
                continue;
            }

            // Default to white
            colors_to_allocate[i].x = 1.f;
            colors_to_allocate[i].y = 1.f;
            colors_to_allocate[i].z = 1.f;
        }

        glBufferData(GL_ARRAY_BUFFER, r->N * sizeof(struct reb_vec3df), colors_to_allocate, GL_STATIC_DRAW);
    }

    // this only does something for WHFAST
    reb_simulation_synchronize(r_copy);

    // Update data on GPU 
    for (size_t i=0;i<N;i++){
        struct reb_particle p = r_copy->particles[i];
        data->particle_data[i].x  = (float)p.x;
        data->particle_data[i].y  = (float)p.y;
        data->particle_data[i].z  = (float)p.z;
        data->particle_data[i].r  = (float)p.r;
    }
    // Only advance breadcrumb index if simulation has advanced
    if (r->steps_done != data->breadcrumb_last_steps_done){
        if (r->steps_done < data->breadcrumb_last_steps_done){
            // Something strange is happening. New simulation?
            reb_display_clear_particle_data(data);
        }
        data->breadcrumb_last_steps_done = r->steps_done;
        data->breadcrumb_current_index = (data->breadcrumb_current_index+1) % data->breadcrumb_N_allocated;
    }

    if (data->s.wire && N>1){
        struct reb_particle com = r_copy->particles[0];
        for (size_t i=1;i<N;i++){
            struct reb_particle p = r_copy->particles[i];
            data->orbit_data[i-1].x  = (float)com.x;
            data->orbit_data[i-1].y  = (float)com.y;
            data->orbit_data[i-1].z  = (float)com.z;
            struct reb_orbit o = reb_orbit_from_particle(r_copy->G, p,com);
            data->orbit_data[i-1].a = (float)o.a;
            data->orbit_data[i-1].e = (float)o.e;
            data->orbit_data[i-1].f = (float)o.f;
            data->orbit_data[i-1].omega = (float)o.omega;
            data->orbit_data[i-1].Omega = (float)o.Omega;
            data->orbit_data[i-1].inc = (float)o.inc;
            com = reb_particle_com_of_pair(p,com);
        }
    }
    if (N>0){
        // Fill memory (but not resize)
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        glBufferSubData(GL_ARRAY_BUFFER, data->breadcrumb_current_index*N*sizeof(struct reb_vec4df), N*sizeof(struct reb_vec4df), data->particle_data);
        if (data->s.spheres==1 || data->s.spheres==2){
            glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer_current);
            glBufferSubData(GL_ARRAY_BUFFER, 0, N*sizeof(struct reb_vec4df), data->particle_data);
        }

        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
        glBufferSubData(GL_ARRAY_BUFFER, data->breadcrumb_current_index*(N-1)*sizeof(struct reb_orbit_opengl), (N-1)*sizeof(struct reb_orbit_opengl), data->orbit_data);
        if (data->s.wire){
            glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer_current);
            glBufferSubData(GL_ARRAY_BUFFER, 0, (N-1)*sizeof(struct reb_orbit_opengl), data->orbit_data);
        }
    }

    // Do actual drawing
    double ratio = (double)width/(double)height;
    glViewport(0,0,width,height);
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT );

#ifndef __EMSCRIPTEN__
    glPointSize(15.*data->retina);
#endif

    // Precalculate matrices 
    struct reb_mat4df projection = reb_mat4df_ortho( -1.6*ratio, 1.6*ratio, -1.6,1.6, -2.5,2.5);
    struct reb_mat4df view = data->s.view;
    if (data->s.reference!=SIZE_MAX){
        struct reb_particle p = data->r_copy->particles[data->s.reference];
        view = reb_mat4df_translate(view, -p.x, -p.y, -p.z);
    }

    if(data->s.drawglow) {
        // Check if the render texture dimensions match, if not, regenerate the textures to the new dimensions
        if(data->reb_post_processing.post_processing_texture_w != width || data->reb_post_processing.post_processing_texture_h != height) {
            data->reb_post_processing.post_processing_texture_w = width;
            data->reb_post_processing.post_processing_texture_h = height;

            glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[0]);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

            glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[1]);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
        }

        // The first frame buffer is where the unfiltered scene will be drawn.
        glBindFramebuffer(GL_FRAMEBUFFER, data->reb_post_processing.post_processing_frame_buffers[0]);
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
        glViewport(0, 0, width, height);
    }

    for (int i=-data->s.ghostboxes*data->r_copy->N_ghost_x;i<=data->s.ghostboxes*data->r_copy->N_ghost_x;i++){
        for (int j=-data->s.ghostboxes*data->r_copy->N_ghost_y;j<=data->s.ghostboxes*data->r_copy->N_ghost_y;j++){
            for (int k=-data->s.ghostboxes*data->r_copy->N_ghost_z;k<=data->s.ghostboxes*data->r_copy->N_ghost_z;k++){
                struct reb_vec6d gb = reb_boundary_get_ghostbox(data->r_copy, i,j,k);
                struct reb_mat4df model = reb_mat4df_translate(reb_mat4df_identity(), gb.x, gb.y, gb.z);
                { // Particles
                    struct reb_mat4df mvp = reb_mat4df_multiply(projection, reb_mat4df_multiply(view, model));
                    if (data->s.wire==2){
                        // Orbit Planes
                        glDisable(GL_CULL_FACE);
                        glUseProgram(data->shader_plane.program);
                        glUniformMatrix4fv(data->shader_plane.mvp_location, 1, GL_TRUE, (GLfloat*) mvp.m);
                        glBindVertexArray(data->shader_plane.particle_vao_current);
                        glUniform1i(data->shader_plane.vertex_count_location, data->shader_plane.vertex_count);
                        reb_glDrawArraysInstanced(GL_TRIANGLES, 0, data->shader_plane.vertex_count, N-1);
                        glBindVertexArray(0);
                        glEnable(GL_CULL_FACE);
                    }

                    if (data->s.wire>=1){
                        if(data->s.drawwiresthick == 1) {
                            // Use the shader program and draw all of the allocated orbits
                            glDisable(GL_CULL_FACE);
                            glUseProgram(data->shader_thick_orbit.shader_program);
                            glUniformMatrix4fv(data->shader_thick_orbit.mvp_location, 1, GL_TRUE, (GLfloat*)mvp.m);
                            glUniformMatrix4fv(data->shader_thick_orbit.v_location, 1, GL_FALSE, (GLfloat*)view.m);
                            glUniformMatrix4fv(data->shader_thick_orbit.p_location, 1, GL_FALSE, (GLfloat*)projection.m);
                            glUniform1i(data->shader_thick_orbit.focused_orbit_section_location, data->shader_thick_orbit.focused_orbit_section_to_use);
                            glUniform1i(data->shader_thick_orbit.anomaly_sample_count_location, data->shader_thick_orbit.anomaly_sample_count_to_use);

                            glBindVertexArray(data->shader_thick_orbit.thick_orbit_vertex_array);
                            reb_glDrawElementsInstanced(GL_TRIANGLES, (data->shader_thick_orbit.anomaly_sample_count_to_use - 1) * 2 * 3, GL_UNSIGNED_INT, NULL, N - 1);
                            glEnable(GL_CULL_FACE);
                        } else {
                            // Orbits
                            glUseProgram(data->shader_orbit.program);
                            glUniformMatrix4fv(data->shader_orbit.mvp_location, 1, GL_TRUE, (GLfloat*) mvp.m);
                            glUniform1i(data->shader_orbit.breadcrumb_N_location, data->breadcrumb_N_allocated);
                            glUniform1i(data->shader_orbit.vertex_count_location, data->shader_orbit.vertex_count);
                            if (data->breadcrumb_N_allocated>1){
                                glBindVertexArray(data->shader_orbit.particle_vao);
                                glUniform1i(data->shader_orbit.N_location, N-1);
                                glUniform1i(data->shader_orbit.current_index_location, data->breadcrumb_current_index);
                                reb_glDrawArraysInstanced(GL_LINE_STRIP, 0, data->shader_orbit.vertex_count, data->breadcrumb_N_allocated*(N-1));
                            }else{
                                glBindVertexArray(data->shader_orbit.particle_vao_current);
                                glUniform1i(data->shader_orbit.N_location, 0);
                                reb_glDrawArraysInstanced(GL_LINE_STRIP, 0, data->shader_orbit.vertex_count, N-1);
                            }
                            glBindVertexArray(0);
                        }
                    }

                    if(data->s.spheres==1||data->s.spheres==2){
                        // Solid Spheres
                        glEnable(GL_DEPTH_TEST);
                        glUseProgram(data->shader_sphere.program);
                        glUniformMatrix4fv(data->shader_sphere.mvp_location, 1, GL_TRUE, (GLfloat*) mvp.m);
                        if (data->breadcrumb_N_allocated>1){
                            glBindVertexArray(data->shader_sphere.particle_vao);
                            reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 800, N*data->breadcrumb_N_allocated);
                        }else{
                            glBindVertexArray(data->shader_sphere.particle_vao_current);
                            reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 800, N);
                        }
                        glBindVertexArray(0);
                        glDisable(GL_DEPTH_TEST);
                    }

                    if(data->s.spheres%2==0){
                        glUseProgram(data->shader_point.program);
                        glBindVertexArray(data->shader_point.particle_vao);
                        glUniformMatrix4fv(data->shader_point.mvp_location, 1, GL_TRUE, (GLfloat*) mvp.m);
                        glUniform1i(data->shader_point.breadcrumb_N_location, data->breadcrumb_N_allocated);
                        if (data->breadcrumb_N_allocated>1){
                            glUniform4f(data->shader_point.color_location, 1.,1.,1.,0.8);
                            glUniform1i(data->shader_point.N_location, N);
                            glUniform1i(data->shader_point.current_index_location, data->breadcrumb_current_index);
                            glDrawArrays(GL_POINTS, 0, N*data->breadcrumb_N_allocated);
                        }
                        // Points
                        glUniform4f(data->shader_point.color_location, 1.,1.,0.,0.8);
                        glUniform1i(data->shader_point.N_location, 0);
                        glDrawArrays(GL_POINTS, N*data->breadcrumb_current_index, N);
                        glBindVertexArray(0);
                    }
                }
                { // Box
                    glUseProgram(data->shader_box.program);
                    struct reb_mat4df boxmodel =  model;
                    if (data->r_copy->boundary == REB_BOUNDARY_NONE){
                        struct reb_vec3df scale = reb_mat4df_get_scale(view); // Extract scale from view matrix so it can be undone
                        glBindVertexArray(data->shader_box.cross_vao);
                        boxmodel = reb_mat4df_scale(boxmodel, 1./scale.x, 1./scale.y, 1./scale.z);
                    }else{
                        glBindVertexArray(data->shader_box.box_vao);
                        boxmodel = reb_mat4df_scale(boxmodel, data->r_copy->root_size*(double)r->N_root_x/2., data->r_copy->root_size*(double)r->N_root_y/2., data->r_copy->root_size*(double)r->N_root_z/2.);
                    }
                    struct reb_mat4df mvp = reb_mat4df_multiply(projection, reb_mat4df_multiply(view, boxmodel));
                    glUniformMatrix4fv(data->shader_box.mvp_location, 1, GL_TRUE, (GLfloat*) mvp.m);
                    glUniform4f(data->shader_box.color_location, 1.,0.,0.,1.);
                    if (data->r_copy->boundary == REB_BOUNDARY_NONE){
                        glDrawArrays(GL_LINES, 0, 2); // x-direction
                        glUniform4f(data->shader_box.color_location, 0.,1.,0.,1.);
                        glDrawArrays(GL_LINES, 2, 2); // y-direction
                        glUniform4f(data->shader_box.color_location, 0.,0.,1.,1.);
                        glDrawArrays(GL_LINES, 4, 2); // z-direction
                        glUniform4f(data->shader_box.color_location, 1.,0.,0.,1.);
                    }else{
                        // x-direction box edges
                        glDrawArrays(GL_LINES, 0, 2);
                        glDrawArrays(GL_LINES, 4, 2);
                        glDrawArrays(GL_LINES, 8, 2);
                        glDrawArrays(GL_LINES, 12, 2);

                        // y-direction box edges
                        glUniform4f(data->shader_box.color_location, 0.,1.,0.,1.);
                        glDrawArrays(GL_LINES, 2, 2);
                        glDrawArrays(GL_LINES, 6, 2);
                        glDrawArrays(GL_LINES, 10, 2);
                        glDrawArrays(GL_LINES, 14, 2);

                        // z-direction box edges
                        glUniform4f(data->shader_box.color_location, 0.,0.,1.,1.);
                        glDrawArrays(GL_LINES, 16, 8);

                        glUniform4f(data->shader_box.color_location, 1.,0.,0.,1.); // back to red
                    }
                    glBindVertexArray(0);
                }
            }
        }
    }

    if(data->s.drawglow) {
        // Bind to next frame buffer to render the horizontal blur to
        glBindFramebuffer(GL_FRAMEBUFFER, data->reb_post_processing.post_processing_frame_buffers[1]);
        glViewport(0,0,width,height);
        glClear(GL_COLOR_BUFFER_BIT);

        // Activate the color texture that was rendered to by the first frame buffer.
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[0]);

        // Apply a horizontal blur and store that in frame buffer #2
        glBindVertexArray(data->reb_post_processing.fullscreen_quad_vao);
        glUseProgram(data->reb_post_processing.horizontal_blur_shader);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Time to draw to the screen now!
        // Dont need to call glClear, already done before binding the frame buffers.
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0,0,width,height);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[0]);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[1]);

        // Draw with the vertical blur and combine shader now.
        glUseProgram(data->reb_post_processing.vertical_blur_and_combine_shader);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    }

    // Ruler
    if (data->s.onscreentext){ 
        glUseProgram(data->shader_box.program);
        glBindVertexArray(data->shader_box.ruler_vao);
        glUniform4f(data->shader_box.color_location, 1.,1.,1.,1.);
        struct reb_vec3df scale3 = reb_mat4df_get_scale(view); // Extract scale from view matrix so it can be undone
        float scaley = powf(10.,floor(log10f(3./scale3.y))); // nearest power of 10, factor of 3. determines wrapping
        if (5.*scaley<3./scale3.y) { scaley*=5;}
        if (2.*scaley<3./scale3.y) { scaley*=2;}

        struct reb_mat4df ruler_mvp =  reb_mat4df_identity();
        ruler_mvp = reb_mat4df_translate(ruler_mvp, 1.-30./width, 0, 0);
        ruler_mvp = reb_mat4df_scale(ruler_mvp, 15.0/width, 0.3125*scale3.y*scaley, 1);  // 0.3125 comes from b and t values in projection matrix
        glUniformMatrix4fv(data->shader_box.mvp_location, 1, GL_TRUE, (GLfloat*) ruler_mvp.m);
        glDrawArrays(GL_LINES, 0, 6);
        glBindVertexArray(0);

        // Text
        char str[256];
        float val[200] = {0.};
        float char_size = data->retina*16.; // px per char
        float scale = 2.*char_size/height; // size of one char in screen coordinates
        glUseProgram(data->shader_simplefont.program);
        glBindVertexArray(data->shader_simplefont.vao);
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(data->shader_simplefont.texture_location, 0);
        glBindTexture(GL_TEXTURE_2D,data->shader_simplefont.texture);
        float screen_aspect = (float)height/(float)width;
        glUniform1f(data->shader_simplefont.screen_aspect_location, screen_aspect);
        glBindBuffer(GL_ARRAY_BUFFER, data->shader_simplefont.charval_buffer);


        // Ruler
        if (scaley >= 1000. || scaley<=0.01){
            sprintf(str, "%.0e", scaley);
        }else{
            sprintf(str, "%.*f", MAX(0,1-(int)log10f(scaley)), scaley);
        }
        float ruler_height = strlen(str)*0.75*scale; 
        glUniform2f(data->shader_simplefont.pos_location, 1.-31./width,-ruler_height/2.);
        glUniform1f(data->shader_simplefont.ypos_location, 0);
        glUniform1i(data->shader_simplefont.rotation_location, 1);
        glUniform1f(data->shader_simplefont.scale_location, scale);
        glUniform1f(data->shader_simplefont.aspect_location, 0.75);
        int j = convertLine(str,val);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

#ifdef __EMSCRIPTEN__
    }
#else // __EMSCRIPTEN__
      // Logo
    char_size = data->retina*4.; // px per char
    scale = 2.*char_size/height; // size of one char in screen coordinates
    float logo_width = 42.0*0.5 *scale*screen_aspect;     //  41=num char, 0.5=aspect
    float logo_height = 26.0 *scale;         //  26=num char
    glUniform2f(data->shader_simplefont.pos_location, -1.,-1.+logo_height);
    glUniform1f(data->shader_simplefont.aspect_location, 0.5);
    glUniform1i(data->shader_simplefont.rotation_location, 0);
    glUniform1f(data->shader_simplefont.scale_location, scale);
    for (size_t i=0;i<sizeof(reb_logo)/sizeof(reb_logo[0]);i++){
        int j = convertLine(reb_logo[i],val);
        glUniform1f(data->shader_simplefont.ypos_location, (float)i);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
    }

    // Status text
    char_size = data->retina*16.; // px per char
    scale = 2.*char_size/height; // size of one char in screen coordinates

    int ypos = 1;
    glUniform2f(data->shader_simplefont.pos_location, -1+logo_width,-1.+logo_height);
    glUniform1f(data->shader_simplefont.aspect_location,0.75);
    glUniform1f(data->shader_simplefont.scale_location, scale);

    glUniform1f(data->shader_simplefont.ypos_location, ypos++);
    sprintf(str,"REBOUND v%s",reb_version_str);
    j = convertLine(str,val);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
    reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

    if (data->r_copy->status == REB_STATUS_RUNNING){
        sprintf(str, "Simulation is running  ");
    }else if (data->r_copy->status == REB_STATUS_PAUSED){
        sprintf(str, "Simulation is paused   ");
    }else if (data->r_copy->status <= REB_STATUS_SINGLE_STEP){
        if (data->r_copy->status == REB_STATUS_SINGLE_STEP){
            sprintf(str, "Integrating 1 step");
        }else{
            sprintf(str, "Integrating %d steps",REB_STATUS_SINGLE_STEP - data->r_copy->status + 1);
        }
    }
    glUniform1f(data->shader_simplefont.ypos_location, ypos++);
    j = convertLine(str,val);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
    reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

    if (!r_copy->display_settings){
        sprintf(str, "Press h for help ");
    }else{
        sprintf(str, "User interaction disabled");
    }
    glUniform1f(data->shader_simplefont.ypos_location, ypos++);
    j = convertLine(str,val);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
    reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);


    sprintf(str, "N = %zu ",data->r_copy->N);
    glUniform1f(data->shader_simplefont.ypos_location, ypos++);
    j = convertLine(str,val);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
    reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

    glUniform1f(data->shader_simplefont.ypos_location, ypos++);
    if (data->r_copy->OMEGA){
        sprintf(str, "t = %f [orb]  ", data->r_copy->t*data->r_copy->OMEGA/2./M_PI);
    }else{
        sprintf(str, "t = %f  ", data->r_copy->t);
    }
    j = convertLine(str,val);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
    reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);

}
glBindVertexArray(0);
glBindTexture(GL_TEXTURE_2D,0);
if (data->s.onscreenhelp){ // On screen help
    glUseProgram(data->shader_simplefont.program);
    glBindVertexArray(data->shader_simplefont.vao);
    glBindTexture(GL_TEXTURE_2D,data->shader_simplefont.texture);
    glUniform2f(data->shader_simplefont.pos_location, -0.67,0.7);
    glUniform1f(data->shader_simplefont.aspect_location, 0.75);
    glUniform1f(data->shader_simplefont.screen_aspect_location, 1./ratio);
    glUniform1i(data->shader_simplefont.rotation_location, 0);
    glUniform1f(data->shader_simplefont.scale_location, 0.035);
    glBindBuffer(GL_ARRAY_BUFFER, data->shader_simplefont.charval_buffer);
    float val[200] = {0.};
    for (size_t i=0;i<sizeof(onscreenhelp)/sizeof(onscreenhelp[0]);i++){
        int j = convertLine(onscreenhelp[i],val);
        glUniform1f(data->shader_simplefont.ypos_location, (float)i);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(val), val);
        reb_glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, j);
    }
}
#endif // __EMSCRIPTEN__

    glfwSwapBuffers(data->window);
    glfwPollEvents();
}

#ifdef __EMSCRIPTEN__
EM_BOOL reb_render_frame_emscripten(double time, void* p){
    struct reb_simulation* r = (struct reb_simulation*)p;
    struct reb_display_data* data = r->display_data;
    if (!data){
        return EM_TRUE;
    }
    if (data->s.pause){
        return EM_TRUE;
    }
    reb_render_frame(data);
    EM_ASM({
            var overlay = document.getElementById("overlay");
            if ($0){
            overlay.style.display = "none";
            }else{
            overlay.style.display = "block";
            }
            }, !data->s.onscreentext);
    if (data->s.onscreentext){ 
        char str[10240] = "\0";
        char line[1024];
        sprintf(line,"<div class=\"reboundlogo\"></div>REBOUND v%s<br />",reb_version_str);
        strlcat(str, line, 10240);
        if (data->connection_status>=0){
            if (data->r_copy->status == REB_STATUS_RUNNING){
                sprintf(line, "Simulation is running<br />");
            }else if (data->r_copy->status == REB_STATUS_PAUSED){
                sprintf(line, "Simulation is paused<br />");
            }else if (data->r_copy->status == REB_STATUS_SCREENSHOT_READY){
                sprintf(line, "Screenshot ready<br />");
            }else if (data->r_copy->status == REB_STATUS_SCREENSHOT){
                sprintf(line, "Taking screenshot<br />");
            }else if (data->r_copy->status == REB_STATUS_SUCCESS){
                sprintf(line, "Simulation ready<br />");
            }else if (data->r_copy->status == REB_STATUS_USER){
                sprintf(line, "Simulation canceled<br />");
            }else if (data->r_copy->status > 0){
                sprintf(line, "Simulation error occurred<br />");
            }else if (data->r_copy->status <= REB_STATUS_SINGLE_STEP){
                if (data->r_copy->status == REB_STATUS_SINGLE_STEP){
                    sprintf(line, "Integrating 1 step<br />");
                }else{
                    sprintf(line, "Integrating %d steps<br />",REB_STATUS_SINGLE_STEP - data->r_copy->status + 1);
                }
            }
            strlcat(str, line, 10240);
            sprintf(line, "N = %zu<br />",data->r_copy->N);
            strlcat(str, line, 10240);
            sprintf(line, "t = %g<br />",data->r_copy->t);
            strlcat(str, line, 10240);
            sprintf(line, "steps/s = %g<br />",1./data->r_copy->walltime_last_steps);
            strlcat(str, line, 10240);
            if (!data->r_copy->display_settings){
                strlcat(str, "Press h or click for help<br />", 10240);
            }else{
                strlcat(str, "User interaction disabled<br />", 10240);
            }
            reb_overlay_update(str, data->r_copy->status);
        }else{
            sprintf(line, "Unable to connect. Server might have shut down.");
            strlcat(str, line, 10240);
            reb_overlay_update(str, 10);
        }
    }
    data->s.onscreenhelp = reb_overlay_help_show(data->s.onscreenhelp);
    if (data->s.onscreenhelp){ 
        char str[10240] = "\0";
        for (int i=0;i<sizeof(onscreenhelp)/sizeof(onscreenhelp[0]);i++){
            strlcat(str, onscreenhelp[i], 10240);
            strlcat(str, "<br />", 10240);
            EM_ASM({
                    var overlaytext = document.getElementById("overlaytext-help");
                    if (overlaytext){
                    overlaytext.innerHTML = UTF8ToString($0);
                    }}, str);
        }
    }
    if (data->take_one_screenshot){
        EM_ASM_PTR({
                var canvas = document.getElementById('canvas');
                var link = document.createElement("a");
                link.download = "screenshot.png";
                link.href = canvas.toDataURL();
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                delete link;
                });
        data->take_one_screenshot = 0;
    }
    if (data->r_copy->status == REB_STATUS_SCREENSHOT && !data->screenshot){
        data->screenshot = EM_ASM_PTR({
                var canvas = document.getElementById('canvas');
                return stringToNewUTF8(canvas.toDataURL());
                });
        if (!data->screenshot){
            printf("Error, screenshot not successful.");
        }
        data->r->status = REB_STATUS_SCREENSHOT_READY; // changing main simulation as r_copy will be overwritten
    }
    return EM_TRUE;
}
#endif


void reb_display_init(struct reb_simulation * const r){
    struct reb_display_data* data = r->display_data;
    if (!glfwInit()){
        reb_simulation_error(r, "GLFW initialization failed.");
        return;
    }

    glfwSetErrorCallback(reb_glfw_error_callback);
#ifdef __EMSCRIPTEN__
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#else
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif

    GLFWwindow*  window = glfwCreateWindow(700, 700, "rebound", NULL, NULL);
    if (!window){
        reb_simulation_error(r,"GLFW window creation failed.");
        return;
    }

    glfwMakeContextCurrent(window);

#ifndef __EMSCRIPTEN__
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
#endif // __EMSCRIPTEN__

    glfwSetWindowUserPointer(window,data); 

    // Default parameters
    reb_display_settings_init(r, &r->display_data->s);
    { // Check if we have a retina display
        int wwidth, wheight, fwidth, fheight;
        glfwGetWindowSize(window, &wwidth, &wheight);
        glfwGetFramebufferSize(window, &fwidth, &fheight);
        data->retina = (double)fwidth/(double)wwidth;
    }
    data->window            = window;
    data->breadcrumb_current_index= 0;
    data->breadcrumb_N_allocated  = 0;

    glfwSetKeyCallback(window,reb_display_keyboard);
    glfwSetScrollCallback(window,reb_display_scroll);
    glfwGetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS);
    glfwSetMouseButtonCallback(window, reb_display_mouse_button);
    glfwSetCursorPosCallback(window, reb_display_cursor);
    glfwSetWindowSizeCallback(window, reb_display_resize);
    glDepthMask(GL_TRUE);
#ifndef __EMSCRIPTEN__
    glEnable(GL_MULTISAMPLE); 
#endif
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendEquation(GL_FUNC_ADD);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT); 

    {
        // SIMPLEFONT shader
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec2 vp;\n"
        "in float charpos;\n"
        "uniform float scale;\n"
        "uniform float aspect;\n"
        "uniform float screen_aspect;\n"
        "uniform float ypos;\n"
        "uniform vec2 pos;\n"
        "uniform int rotation;\n"
        "in vec2 charval;\n"
        "in vec2 texcoord;\n"
        "out vec2 Texcoord;\n"
        "void main() {\n"
        "  if (rotation==0) {\n"
        "    gl_Position = vec4(pos.x+screen_aspect*scale*(vp.x+charpos*aspect),pos.y+scale*(vp.y-ypos),0.,1.);\n"
        "  }else{\n"
        "    gl_Position = vec4(pos.x+screen_aspect*scale*(-vp.y),pos.y+scale*(vp.x-ypos+charpos*aspect),0.,1.);\n"
        "  }\n"
        "  Texcoord = vec2((charval.s+texcoord.s)/16.,(charval.t+texcoord.t)/16.00);\n"
        "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "out vec4 outcolor;\n"
        "uniform sampler2D tex;\n"
        "in vec2 Texcoord;\n"
        "void main() {\n"
        "  outcolor =   vec4(0.5,0.5,0.5,(1.-texture(tex, Texcoord).r)); \n"
        "}\n";

        data->shader_simplefont.program = loadShader(vertex_shader, fragment_shader);
        data->shader_simplefont.ypos_location = glGetUniformLocation(data->shader_simplefont.program, "ypos");
        data->shader_simplefont.pos_location = glGetUniformLocation(data->shader_simplefont.program, "pos");
        data->shader_simplefont.scale_location = glGetUniformLocation(data->shader_simplefont.program, "scale");
        data->shader_simplefont.rotation_location = glGetUniformLocation(data->shader_simplefont.program, "rotation");
        data->shader_simplefont.texture_location = glGetUniformLocation(data->shader_simplefont.program, "tex");
        data->shader_simplefont.aspect_location = glGetUniformLocation(data->shader_simplefont.program, "aspect");
        data->shader_simplefont.screen_aspect_location = glGetUniformLocation(data->shader_simplefont.program, "screen_aspect");

        glUseProgram(data->shader_simplefont.program);
        glGenVertexArrays(1, &data->shader_simplefont.vao);
        glBindVertexArray(data->shader_simplefont.vao);
        GLuint sfvp = glGetAttribLocation(data->shader_simplefont.program,"vp");
        glEnableVertexAttribArray(sfvp);
        GLuint sftexcoordp = glGetAttribLocation(data->shader_simplefont.program,"texcoord");
        glEnableVertexAttribArray(sftexcoordp);
        GLuint charval_location = glGetAttribLocation(data->shader_simplefont.program,"charval");
        glEnableVertexAttribArray(charval_location);
        GLuint charpos_location = glGetAttribLocation(data->shader_simplefont.program,"charpos");
        glEnableVertexAttribArray(charpos_location);
        float simplefont_data[] = {
            0., 0., 0., 1., 
            0., 1., 0., 0.0,
            1., 0., 1., 1.,
            1., 1., 1., 0.0
        };
        GLuint simplefont_buffer;
        glGenBuffers(1, &simplefont_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, simplefont_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(simplefont_data), simplefont_data, GL_STATIC_DRAW);

        glVertexAttribPointer(sfvp, 2, GL_FLOAT, GL_FALSE, sizeof(float)*4, NULL);
        glVertexAttribPointer(sftexcoordp, 2, GL_FLOAT, GL_FALSE, sizeof(float)*4, (void *)(sizeof(float)*2));

        GLuint charpos_buffer;
        glGenBuffers(1, &charpos_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, charpos_buffer);
        float charpos[100];
        for(int i=0;i<100;i++){
            charpos[i] = (float)i;
        }
        glBufferData(GL_ARRAY_BUFFER, sizeof(charpos), charpos, GL_STATIC_DRAW);
        glVertexAttribPointer(charpos_location, 1, GL_FLOAT, GL_FALSE, sizeof(float)*1, NULL);

        glGenBuffers(1, &data->shader_simplefont.charval_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, data->shader_simplefont.charval_buffer);
        glBufferData(GL_ARRAY_BUFFER, 200*sizeof(float), NULL,GL_DYNAMIC_DRAW);
        glVertexAttribPointer(charval_location, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, NULL);

        reb_glVertexAttribDivisor(sfvp, 0); 
        reb_glVertexAttribDivisor(sftexcoordp, 0); 
        reb_glVertexAttribDivisor(charpos_location,1);
        reb_glVertexAttribDivisor(charval_location,1);

        glBindVertexArray(0);

        // Load simplefont
        unsigned char image[65536];
        for(int i=0;i<8192;i++){
            unsigned char byte = simplefont[i];
            for(int b=0;b<8;b++){
                image[i*8+b] = ((byte >> b) & 0x01)? 0xFF : 0x00;
            }
        }
        glGenTextures(1, &data->shader_simplefont.texture);
        glBindTexture(GL_TEXTURE_2D,data->shader_simplefont.texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, 256, 256, 0, GL_RED, GL_UNSIGNED_BYTE, image);
        glBindTexture(GL_TEXTURE_2D,0);

    }

    {
        // POINT shader
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec3 vp;\n"
        "uniform mat4 mvp;\n"
        "uniform vec4 vc;\n"
        "out vec4 color;\n"
        "uniform int current_index;\n"
        "uniform int N;\n"
        "uniform int breadcrumb_N;\n"
        "void main() {\n"
        "  gl_Position = mvp*vec4(vp, 1.0);\n"
        "  gl_Position.z = 0.;\n" // no clipping
        "  gl_PointSize = 15.0f;\n"
        "  color = vc;\n"
        "  float age = float( (gl_VertexID/N - current_index + 2*breadcrumb_N -1)%breadcrumb_N +1 )/float(breadcrumb_N);\n"
        "  if (N == 0) age = 1.0;\n"
        "  color = vec4(vc.xyz,age*vc.a);\n"
        "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "in vec4 color;\n"
        "out vec4 outcolor;\n"
        "void main() {\n"
        "  vec2 rel = gl_PointCoord.st;\n"
        "  rel.s -=0.5f;\n"
        "  rel.t -=0.5f;\n"
        "  if (length(rel)>0.25f){\n"
        "     outcolor = vec4(0.f,0.f,0.f,0.f); \n"
        "  }else{\n"
        "     vec4 cmod = color;\n"
        "     cmod.a*= min(1.,1.-4.*(length(rel)/0.25-0.75));\n"
        "     outcolor = cmod;\n"
        "  }\n"
        "}\n";

        data->shader_point.program = loadShader(vertex_shader, fragment_shader);
        data->shader_point.mvp_location = glGetUniformLocation(data->shader_point.program, "mvp");
        data->shader_point.color_location = glGetUniformLocation(data->shader_point.program, "vc");
        data->shader_point.current_index_location = glGetUniformLocation(data->shader_point.program, "current_index");
        data->shader_point.breadcrumb_N_location = glGetUniformLocation(data->shader_point.program, "breadcrumb_N");
        data->shader_point.N_location = glGetUniformLocation(data->shader_point.program, "N");

        glUseProgram(data->shader_point.program);
        glGenVertexArrays(1, &data->shader_point.particle_vao);
        glBindVertexArray(data->shader_point.particle_vao);
        GLuint pvp = glGetAttribLocation(data->shader_point.program,"vp");
        glEnableVertexAttribArray(pvp);

        glGenBuffers(1, &data->particle_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
        glVertexAttribPointer(pvp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*4, NULL);
        glBindVertexArray(0);
    }

    {
        // BOX shader
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec3 vp;\n"
        "uniform mat4 mvp;\n"
        "uniform vec4 vc;\n"
        "out vec4 color;\n"
        "void main() {\n"
        "  gl_Position = mvp*vec4(vp, 1.0);\n"
        "  gl_Position.z = 0.;\n" // no clipping
        "  color = vc;\n"
        "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "in vec4 color;\n"
        "out vec4 outcolor;\n"
        "void main() {\n"
        "  outcolor = color;\n"
        "}\n";

        data->shader_box.program = loadShader(vertex_shader, fragment_shader);
        data->shader_box.mvp_location = glGetUniformLocation(data->shader_box.program, "mvp");
        data->shader_box.color_location = glGetUniformLocation(data->shader_box.program, "vc");

        // Create cross mesh
        glUseProgram(data->shader_box.program);
        glGenVertexArrays(1, &data->shader_box.cross_vao);
        glBindVertexArray(data->shader_box.cross_vao);
        GLuint vp = glGetAttribLocation(data->shader_box.program,"vp");
        glEnableVertexAttribArray(vp);

        float cross_data[18] = {
            -0.04,0.,0., +0.04,0.,0., 0.,-0.04,0., 0.,+0.04,0., 0.,0.,-0.04, 0.,0.,+0.04,
        };
        GLuint cross_buffer;
        glGenBuffers(1, &cross_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, cross_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(cross_data), cross_data, GL_STATIC_DRAW);

        glVertexAttribPointer(vp, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glBindVertexArray(0);

        // Create box mesh
        glGenVertexArrays(1, &data->shader_box.box_vao);
        glBindVertexArray(data->shader_box.box_vao);
        glEnableVertexAttribArray(vp);

        float box_data[] = {
            -1,-1,-1, 1,-1,-1, 1,-1,-1, 1, 1,-1, 1, 1,-1, -1, 1,-1, -1, 1,-1, -1,-1,-1,
            -1,-1, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1,-1, 1,
            -1,-1,-1, -1,-1, 1, -1, 1,-1, -1, 1, 1, 1, 1,-1, 1, 1, 1, 1,-1,-1, 1,-1, 1,
        };
        GLuint box_buffer;
        glGenBuffers(1, &box_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, box_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(box_data), box_data, GL_STATIC_DRAW);

        glVertexAttribPointer(vp, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glBindVertexArray(0);

        // Create ruler mesh
        glGenVertexArrays(1, &data->shader_box.ruler_vao);
        glBindVertexArray(data->shader_box.ruler_vao);
        glEnableVertexAttribArray(vp);

        float ruler_data[18] = {
            0.0, -1.0, 0., 
            0.0, 1.0, 0., 
            1.0, 1.0, 0., 
            -1.0, 1.0, 0., 
            1.0, -1.0, 0., 
            -1.0, -1.0, 0., 
        };
        GLuint ruler_buffer;
        glGenBuffers(1, &ruler_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, ruler_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ruler_data), ruler_data, GL_STATIC_DRAW);

        glVertexAttribPointer(vp, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        glBindVertexArray(0);

    }

    {
        // SPHERE shader
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec3 vp;\n"
        "in float sr;\n"
        "in vec3 sp;\n"
        "out vec3 normal;\n"
        "uniform mat4 mvp;\n"
        "void main() {\n"
        "  gl_Position = mvp*(vec4(sr*vp, 1.0)+vec4(sp,0.));\n"
        "  normal = vp;\n"
        "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "out vec4 outcolor;\n"
        "in vec3 normal;\n"
        "void main() {\n"
        "  vec3 lightdir = vec3(1.,1.,1.);\n"
        "  float intensity = 0.5+max(0.,0.5*dot(normalize(lightdir),normalize(normal)));\n"
        "  outcolor = vec4(intensity,intensity,intensity,1.);\n"
        "}\n";

        data->shader_sphere.program = loadShader(vertex_shader, fragment_shader);
        data->shader_sphere.mvp_location = glGetUniformLocation(data->shader_sphere.program, "mvp");

        // Sphere data
        float* sphere_data = malloc(sizeof(float)*3*800);
        int count = 0; // will be 800 by end
        const int ni = 20;
        const int nj = 20;
        for(int i=0;i<ni;i++){
            for(int j=0;j<nj;j++){
                sphere_data[count*3+0] = sinf((float)i/(float)ni*M_PI)*sinf((float)j/(float)nj*2.*M_PI);
                sphere_data[count*3+1] = sinf((float)i/(float)ni*M_PI)*cosf((float)j/(float)nj*2.*M_PI);
                sphere_data[count*3+2] = cosf((float)i/(float)ni*M_PI);
                count++;
                sphere_data[count*3+0] = sinf((float)(i+1)/(float)ni*M_PI)*sinf((float)j/(float)nj*2.*M_PI);
                sphere_data[count*3+1] = sinf((float)(i+1)/(float)ni*M_PI)*cosf((float)j/(float)nj*2.*M_PI);
                sphere_data[count*3+2] = cosf((float)(i+1)/(float)ni*M_PI);
                count++;
            }
        }
        GLuint sphere_vertex_buffer;
        glGenBuffers(1, &sphere_vertex_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, sphere_vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*count, sphere_data, GL_STATIC_DRAW);
        free(sphere_data);


        // Create sphere vao (one for current, one for past)
        glUseProgram(data->shader_sphere.program);
        GLuint svp = glGetAttribLocation(data->shader_sphere.program,"vp");
        GLuint ssp = glGetAttribLocation(data->shader_sphere.program,"sp");
        GLuint ssr = glGetAttribLocation(data->shader_sphere.program,"sr");


        { // current
            glGenVertexArrays(1, &data->shader_sphere.particle_vao_current);
            glBindVertexArray(data->shader_sphere.particle_vao_current);
            glBindBuffer(GL_ARRAY_BUFFER, sphere_vertex_buffer);
            glEnableVertexAttribArray(svp);

            glVertexAttribPointer(svp, 3, GL_FLOAT, GL_FALSE, 0, NULL);

            glEnableVertexAttribArray(ssp);
            glEnableVertexAttribArray(ssr);
            glGenBuffers(1, &data->particle_buffer_current);
            glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer_current);
            glVertexAttribPointer(ssp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*4, NULL);
            glVertexAttribPointer(ssr, 1, GL_FLOAT, GL_FALSE, sizeof(float)*4, (void*)(sizeof(float)*3));

            reb_glVertexAttribDivisor(svp, 0); 
            reb_glVertexAttribDivisor(data->shader_sphere.mvp_location, 0); 
            reb_glVertexAttribDivisor(ssp, 1);
            reb_glVertexAttribDivisor(ssr, 1);
        }
        { // current
            glGenVertexArrays(1, &data->shader_sphere.particle_vao);
            glBindVertexArray(data->shader_sphere.particle_vao);
            glBindBuffer(GL_ARRAY_BUFFER, sphere_vertex_buffer);
            glEnableVertexAttribArray(svp);

            glVertexAttribPointer(svp, 3, GL_FLOAT, GL_FALSE, 0, NULL);

            glEnableVertexAttribArray(ssp);
            glEnableVertexAttribArray(ssr);
            // glGenBuffers(1, &data->particle_buffer); // generated earlier
            glBindBuffer(GL_ARRAY_BUFFER, data->particle_buffer);
            glVertexAttribPointer(ssp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*4, NULL);
            glVertexAttribPointer(ssr, 1, GL_FLOAT, GL_FALSE, sizeof(float)*4, (void*)(sizeof(float)*3));

            reb_glVertexAttribDivisor(svp, 0); 
            reb_glVertexAttribDivisor(data->shader_sphere.mvp_location, 0); 
            reb_glVertexAttribDivisor(ssp, 1);
            reb_glVertexAttribDivisor(ssr, 1);
        }

        glBindVertexArray(0);
    }

    {
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec3 focus;\n"
        "in vec3 aef;\n"
        "in vec3 omegaOmegainc;\n"
        "uniform int current_index;\n"
        "uniform int N;\n"
        "uniform int vertex_count;\n"
        "uniform int breadcrumb_N;\n"
        "out float lin;\n"
        "uniform mat4 mvp;\n"
        "const float M_PI = 3.14159265359;\n"

        "float true_to_eccentric(float f, float e) {\n"
        "   float E = 2.0 * atan(sqrt(1.0 - e) * sin(f / 2.0), sqrt(1.0 + e) * cos(f / 2.0));\n"
        "   if (E < 0.0) E += 2.0 * M_PI;\n"
        "   return E;\n"
        "}\n"

        "void main() {\n"
        "   float a = aef.x;\n"
        "   float e = aef.y;\n"
        "   lin = float(gl_VertexID)/float(vertex_count-1);\n"
        "   float f = aef.z;\n"
        "   float E = true_to_eccentric(f, e) + lin * 2.f * M_PI;\n"
        "   if (e>1.){\n"
        "       float theta_max = acos(-1./e);\n"
        "       f = 0.0001-theta_max+1.9998*lin*theta_max;\n"
        "       lin = sqrt(min(0.5,lin));\n"
        "   }\n"
        "   if (N != 0) {\n"
        "       lin = float( (gl_InstanceID/N - current_index + 2*breadcrumb_N -1)%breadcrumb_N )/float(breadcrumb_N);\n"
        "   }\n"

        "   float omega = omegaOmegainc.x;\n"
        "   float Omega = omegaOmegainc.y;\n"
        "   float inc = omegaOmegainc.z;\n"
        "   float cO = cos(Omega);\n"
        "   float sO = sin(Omega);\n"
        "   float co = cos(omega);\n"
        "   float so = sin(omega);\n"
        "   float ci = cos(inc);\n"
        "   float si = sin(inc);\n"

        // If hyperbolic, use old true anomaly sampling
        "   if(e > 1.f) {\n"
        "       float cf = cos(f);\n"
        "       float sf = sin(f);\n"
        "       float r = a*(1.-e*e)/(1. + e*cos(f));\n"
        "       vec3 pos = vec3(r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci),r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci),+ r*(so*cf+co*sf)*si);\n"
        "       gl_Position = mvp*(vec4(focus+pos, 1.0));\n"
        "       gl_Position.z = 0.;\n" // no clipping
        "       return;\n"
        "   }\n"

        // The formula for getting the position using eccentric anomaly given orbital parameters
        "   float cE = cos(E);\n"
        "   float sE = sin(E);\n"
        "   float x = a*(cE - e);\n"
        "   float y = a*sqrt(1. - e*e)*sE;\n"
        "   vec3 pos = vec3(\n"
        "       cO*(co*x - so*y) - sO*(so*x*ci + co*y*ci),\n"
        "       sO*(co*x - so*y) + cO*(so*x*ci + co*y*ci),\n"
        "       si*(so*x + co*y)\n"
        "   );\n"
        "   gl_Position = mvp*(vec4(focus+pos, 1.0));\n"
        "   gl_Position.z = 0.;\n" // no clipping
        "}\n";

        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "out vec4 outcolor;\n"
        "in float lin;\n"
        "void main() {\n"
        "  outcolor = vec4(1.,1.,1.,sqrt(lin));\n"
        "}\n";

        data->shader_orbit.program = loadShader(vertex_shader, fragment_shader);
        data->shader_orbit.mvp_location = glGetUniformLocation(data->shader_orbit.program, "mvp");
        data->shader_orbit.current_index_location = glGetUniformLocation(data->shader_orbit.program, "current_index");
        data->shader_orbit.vertex_count_location = glGetUniformLocation(data->shader_orbit.program, "vertex_count");
        data->shader_orbit.breadcrumb_N_location = glGetUniformLocation(data->shader_orbit.program, "breadcrumb_N");
        data->shader_orbit.N_location = glGetUniformLocation(data->shader_orbit.program, "N");

        // Original vertex_count (500)
        data->shader_orbit.vertex_count = 500; // higher number = smoother orbits

        // Generate two orbit vao
        glUseProgram(data->shader_orbit.program);
        GLuint ofocusp = glGetAttribLocation(data->shader_orbit.program,"focus");
        GLuint oaefp = glGetAttribLocation(data->shader_orbit.program,"aef");
        GLuint oomegaOmegaincp = glGetAttribLocation(data->shader_orbit.program,"omegaOmegainc");

        { // Current
            glGenVertexArrays(1, &data->shader_orbit.particle_vao_current);
            glBindVertexArray(data->shader_orbit.particle_vao_current);
            glEnableVertexAttribArray(ofocusp);
            glEnableVertexAttribArray(oaefp);
            glEnableVertexAttribArray(oomegaOmegaincp);

            glGenBuffers(1, &data->orbit_buffer_current);
            glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer_current);
            glVertexAttribPointer(ofocusp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, NULL);
            glVertexAttribPointer(oaefp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*3));
            glVertexAttribPointer(oomegaOmegaincp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*6));

            reb_glVertexAttribDivisor(data->shader_orbit.mvp_location, 0); 
            reb_glVertexAttribDivisor(ofocusp, 1);
            reb_glVertexAttribDivisor(oaefp, 1);
            reb_glVertexAttribDivisor(oomegaOmegaincp, 1);
        }

        { // Past
            glGenVertexArrays(1, &data->shader_orbit.particle_vao);
            glBindVertexArray(data->shader_orbit.particle_vao);
            glEnableVertexAttribArray(ofocusp);
            glEnableVertexAttribArray(oaefp);
            glEnableVertexAttribArray(oomegaOmegaincp);

            glGenBuffers(1, &data->orbit_buffer);
            glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer);
            glVertexAttribPointer(ofocusp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, NULL);
            glVertexAttribPointer(oaefp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*3));
            glVertexAttribPointer(oomegaOmegaincp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*6));

            reb_glVertexAttribDivisor(data->shader_orbit.mvp_location, 0); 
            reb_glVertexAttribDivisor(ofocusp, 1);
            reb_glVertexAttribDivisor(oaefp, 1);
            reb_glVertexAttribDivisor(oomegaOmegaincp, 1);
        }

        glBindVertexArray(0);
    }

    { // Thickened Orbit Rendering
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "const float distribution_factor = 0.5;"
        "const float distribution_interp_power = 5.0;\n"

        "const float width_perturbation_factor = .025f;\n"
        "const float M_PI = 3.14159265359;\n"

        "in vec3 trail_color;\n"
        "in vec3 focus;\n"
        "in vec3 aef;\n"
        "in vec3 omegaOmegainc;\n"

        "uniform int anomaly_sample_count;\n"
        "uniform int orbital_section_focused;\n"
        "out float lin;\n"
        "out float alpha_to_use;\n"
        "out float tangentlin;\n"
        "out vec3 trail_color_to_use;\n"
        "uniform mat4 v;\n"
        "uniform mat4 p;\n"
        "uniform mat4 mvp;\n"

        "vec3 rotate_to_orbit_plane(float E, float x, float y) {\n"
        "   float omega = omegaOmegainc.x;\n"
        "   float Omega = omegaOmegainc.y;\n"
        "   float inc = omegaOmegainc.z;\n"
        "   float cE = cos(E);\n"
        "   float sE = sin(E);\n"
        "   float cO = cos(Omega);\n"
        "   float sO = sin(Omega);\n"
        "   float co = cos(omega);\n"
        "   float so = sin(omega);\n"
        "   float ci = cos(inc);\n"
        "   float si = sin(inc);\n"
        "   return vec3(\n"
        "       cO*(co*x - so*y) - sO*(so*x*ci + co*y*ci),\n"
        "       sO*(co*x - so*y) + cO*(so*x*ci + co*y*ci),\n"
        "       si*(so*x + co*y)\n"
        "   );\n"
        "}\n"

        // Converts true anomaly to eccentric anomaly
        "float true_to_eccentric(float f, float e) {\n"
        "   float E = 2.0 * atan(sqrt(1.0 - e) * sin(f / 2.0), sqrt(1.0 + e) * cos(f / 2.0));\n"
        "   if (E < 0.0) E += 2.0 * M_PI;\n"
        "   return E;\n"
        "}\n"

        "void calculate_perturbation(in float even_E, in float E_not, in float sample_id, out float perturbed_E) {\n"
        "   if(even_E == 0.f) { perturbed_E = 0.f; return; }\n"

        // Initial parameters
        "   float E_per_sample = 2.f * M_PI / float(anomaly_sample_count - 1);\n"
        "   float max_multiplier = (E_not + E_per_sample * (sample_id + 1.f)) / even_E;\n"
        "   float min_multiplier = (E_not + E_per_sample * (sample_id - 1.f)) / even_E;\n"
        "   float k_0 = 0.5f;\n" // Determines maximum shift to endpoint
        "   float k_1 = 2.f;\n" // Determines effect falloff rate with distance from apsis.
        "   float shift_direction = 0.f;\n"

        // Distance ratio to PI and 0, 2PI endpoints
        "   float dist_to_endpoint = 0.f;\n"
        "   float norm_even_E = mod(even_E, 2.f * M_PI);\n"
        "   if(norm_even_E < 0.5f * M_PI) { shift_direction = -1.f; dist_to_endpoint = norm_even_E; }\n"
        "   if(norm_even_E >= 0.5f * M_PI) { shift_direction = 1.f; dist_to_endpoint = M_PI - norm_even_E; }\n"
        "   if(norm_even_E >= M_PI) { shift_direction = -1.f; dist_to_endpoint = norm_even_E - M_PI; }\n"
        "   if(norm_even_E >= 1.5f * M_PI) { shift_direction = 1.f; dist_to_endpoint = 2.f * M_PI - norm_even_E; }\n"
        "   float endpoint_dist_ratio = 2.f * dist_to_endpoint / M_PI;\n"
        "   if(endpoint_dist_ratio < 0.01f) { endpoint_dist_ratio = 0.f; }\n"

        // Distance ratio to initial value of eccentric anomaly
        // Shifting should not occur when at initial value because that is where the particle is rendered.
        // Shifting at where the particle is rendered will have the trail be out of sync with the particle.
        "   float dist_to_E_not = 0.f;\n"
        "   if(even_E >= E_not + M_PI) { dist_to_E_not = E_not + 2.f * M_PI - even_E; }\n"
        "   if(even_E < E_not + M_PI) { dist_to_E_not = even_E - E_not; }\n"
        "   float origin_dist_ratio = dist_to_E_not / M_PI;\n"
        "   float effect_multiplier = 1.f;\n"
        "   if(origin_dist_ratio < 0.05f) { effect_multiplier = 0.f; }\n"
        "   if(origin_dist_ratio > 0.05f) { effect_multiplier = 1.f; }\n"

        // The multiplier to interpolate towards based on distance and shifting direction
        "   float multiplier_to_use = 1.f;\n"
        "   if(shift_direction == 1.f) { multiplier_to_use = max_multiplier; }\n"
        "   if(shift_direction == -1.f) { multiplier_to_use = min_multiplier; }\n"

        // If very close to endpoint, dont shift as this can cause overshooting/undershooting creating kinks
        "   if(endpoint_dist_ratio == 0.f) { multiplier_to_use = 1.f; }\n"
        "   float multiplier = 1.f + effect_multiplier * pow(k_0 * (1.f - endpoint_dist_ratio), k_1) * (multiplier_to_use - 1.f);\n"
        "   perturbed_E = multiplier * even_E;\n"
        "}\n"

        // This is a custom implementation of the scaling that the cross product gives.
        // This is used to shrink the endpoint artifacts that the fast 180 degree rotation can cause.
        "float calculate_normal_scale(float x, float k) {\n"
        "   float root_k = sqrt(k);\n"
        "   if(x < 1.f / root_k) { return 1.f - k * pow(x - 1.f / root_k, 2.f); }\n"
        "   if(x >= 1.f / root_k) { return 1.f; }\n"
        "   return 0.f;\n"
        "}\n"

        // If the given vertex is an odd vertex, it will be shifted outward (away from ellipse).
        // If the given vertex is an even vertex, it will be shifted inward (towards ellipse).
        "void get_width_perturbation(in int vertex_index, out float width_perturbation) {\n"
        "   int remainder = vertex_index % 2;\n"
        "   if(remainder == 1) { width_perturbation = width_perturbation_factor; }\n"
        "   if(remainder == 0) { width_perturbation = -width_perturbation_factor; }\n"
        "}\n"

        "void main() {\n"
        "   float a = aef.x;\n"
        "   float e = aef.y;\n"
        "   float f = aef.z;\n"
        "   float lin = 0.f;\n"

        "   int vertex_id_to_sample = gl_VertexID;\n"
        "   if (vertex_id_to_sample % 2 != 0) { vertex_id_to_sample -= 1; }"
        "   lin = float(vertex_id_to_sample) / float(2 * (anomaly_sample_count - 1));\n"
        "   trail_color_to_use = trail_color;\n"

        // Get the initial true anomaly
        // Convert it to eccentric anomaly for the algorithm (due to required sampling symmetry over entire orbit)
        "   float initial_E = true_to_eccentric(f, e);\n"
        "   float even_distributed_E = initial_E + lin * M_PI * 2.f;\n"
        "   float perturbed_E = even_distributed_E;\n"
        
        // Calculate the shifted eccentric anomaly to increase detail at sharp turns and save the result
        // Also calculate the alpha to use here for the fragment shader
        "   if(e > 0.5) { calculate_perturbation(even_distributed_E, initial_E, float(vertex_id_to_sample / 2), perturbed_E); }"
        "   float relative_lin = mod(perturbed_E, 2.f * M_PI) / (2.f * M_PI);\n"
        "   if(orbital_section_focused == 0) { alpha_to_use = sqrt(lin); }\n"
        "   if(orbital_section_focused == 1) { alpha_to_use = (lin * lin) * (0.5f + 0.5f * cos(2.f * M_PI * clamp(relative_lin, 0.f, 1.f))); }\n"
        "   if(orbital_section_focused == 2) { alpha_to_use = (lin * lin) * (0.5f - 0.5f * cos(2.f * M_PI * clamp(relative_lin, 0.f, 1.f))); }\n"

        // The formula for getting the position using eccentric anomaly given orbital parameters.
        "   float radial_perturbation = 0.f;\n"
        "   get_width_perturbation(gl_VertexID, radial_perturbation);\n"
        "   float x = a*(cos(perturbed_E) - e);\n"
        "   float y = a*sqrt(1. - e*e)*sin(perturbed_E);\n"

        // This is the 2D rotation matrix with theta=pi/2 applied onto x, y
        "   float tx = -a*sin(perturbed_E);\n"
        "   float ty = a*sqrt(1. - e*e)*cos(perturbed_E);\n"

        // Obtain the tangent and position on the orbital plane, not the 2D plane
        // Also obtain the "forward" of the plane. This will be used to adjust the endpoint scaling
        "   float tm = sqrt(tx * tx + ty * ty);\n"
        "   vec3 pos = rotate_to_orbit_plane(perturbed_E, x, y);\n"
        "   vec3 plane_forward = rotate_to_orbit_plane(0.f, 1.f, 0.f);\n"
        "   vec3 plane_tangent = rotate_to_orbit_plane(perturbed_E, tx / tm, ty / tm);\n"

        // Cross the forward vector of the view matrix and the orbit's tangent on the orbital plane
        "   vec3 view_forward = normalize(vec3(v[2][0], v[2][1], v[2][2]));\n"
        "   vec3 trail_normal = cross(view_forward, plane_tangent);\n"

        // Get the linearized length from [0, 1] of the trail_normal to apply custom scaling
        // This is an attempt at keeping a uniform width throughout entire trail, but shrinking at ends to create nubs
        // This is to hide the fast 180 degree rotation of the trail that causes ugly artifacts at tangential views of the trail
        "   float trail_normal_length_lin = asin(length(trail_normal)) / (M_PI / 2.f);\n"

        // For the trail normal scaling function, I interpolate between a small effect and large effect.
        // This is to account for the rate of observed change of the ellipse tangent if looking down the "short side" or "long side"
        "   float min_k = 1.f;\n" // This value will be used when looking down "long side" of ellipse
        "   float max_k = 100.f;\n" // This value will be used when looking down "short side" of ellipse
        "   float k_value_x = abs(dot(plane_forward, view_forward));\n"
        "   if(e < 0.5) { k_value_x = 0.5f; }"

        // If the eccentricity is small, the effect is applied at 25%. Otherwise, the effect strength scales with e.
        "   float k_to_use = min_k + max(0.25, e) * pow(k_value_x, 5.f) * (max_k - min_k);\n"
        "   float normal_scale_to_apply = calculate_normal_scale(trail_normal_length_lin, k_to_use);\n"

        // To prevent divisions by zero, check if the length of the trail is above some small value
        // If so, we can apply the normal normally.
        // If not, the default normalized trail normal is just zero.
        "   vec3 trail_normal_to_apply = vec3(0.f, 0.f, 0.f);\n"
        "   if(length(trail_normal) > 0.00001f) {\n"
        "       trail_normal_to_apply = trail_normal / length(trail_normal);\n"
        "       trail_normal_to_apply *= normal_scale_to_apply;\n"
        "   }\n"
        
        // Finally apply the normal and radial perturbations calculated earlier to build the trail.
        // In this code, zooming works by uniform scaling of the view matrix
        // I can just take the length of the view's forward vector, and use it to scale the trail.
        "   float view_scaling = length(vec3(v[2][0], v[2][1], v[2][2]));\n"
        "   float scaled_radial_perturbation = (1.f / view_scaling) * radial_perturbation;\n"

        "   vec3 perturbed_pos = pos + trail_normal_to_apply * (scaled_radial_perturbation / 2.f);\n"
        "   gl_Position = mvp * (vec4(focus+perturbed_pos, 1.0));\n"
        "   gl_Position.z = 0.;\n" // no clipping
        "}\n";

        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;\n"
        "out vec4 outcolor;\n"

        "in vec3 trail_color_to_use;\n"
        "in float alpha_to_use;\n"
        
        "void main() {\n"
        "   outcolor = vec4(trail_color_to_use, alpha_to_use);\n"
        "}\n";
        
        // Rendering Settings
        data->shader_thick_orbit.anomaly_sample_count_to_use = 250;
        data->shader_thick_orbit.focused_orbit_section_to_use = 1;

        // Load the above shaders
        data->shader_thick_orbit.shader_program = loadShader(vertex_shader, fragment_shader);

        // Store the uniform location for the anomaly sample count.
        // The shader will sample the eccentric anomaly, and convert that to a position given the orbital parameters.
        // The anomaly sample count will dictate how many points will be uniformly sampled on the orbit.
        // The total number of vertices will be 2 * anomaly_sample_count as one sample corresponds to an inner and outer vertex.
        glUseProgram(data->shader_thick_orbit.shader_program);
        data->shader_thick_orbit.anomaly_sample_count_location = glGetUniformLocation(data->shader_thick_orbit.shader_program, "anomaly_sample_count");
        data->shader_thick_orbit.focused_orbit_section_location = glGetUniformLocation(data->shader_thick_orbit.shader_program, "orbital_section_focused");
        data->shader_thick_orbit.mvp_location = glGetUniformLocation(data->shader_thick_orbit.shader_program, "mvp");
        data->shader_thick_orbit.v_location = glGetUniformLocation(data->shader_thick_orbit.shader_program, "v");
        data->shader_thick_orbit.p_location = glGetUniformLocation(data->shader_thick_orbit.shader_program, "p");
        
        // Get all vertex attributes that the shader will compute with.
        GLuint trail_colors = glGetAttribLocation(data->shader_thick_orbit.shader_program, "trail_color");
        GLuint ofocusp = glGetAttribLocation(data->shader_thick_orbit.shader_program,"focus");
        GLuint oaefp = glGetAttribLocation(data->shader_thick_orbit.shader_program,"aef");
        GLuint oomegaOmegaincp = glGetAttribLocation(data->shader_thick_orbit.shader_program,"omegaOmegainc");
        
        // Create the vertex array and enable each of the attributes retrieved above.
        glGenVertexArrays(1, &data->shader_thick_orbit.thick_orbit_vertex_array);
        glBindVertexArray(data->shader_thick_orbit.thick_orbit_vertex_array);
        glEnableVertexAttribArray(ofocusp);
        glEnableVertexAttribArray(oaefp);
        glEnableVertexAttribArray(oomegaOmegaincp);
        
        // This Vertex Array Object will point to the "present" orbit buffer for drawing.
        glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer_current);
        glVertexAttribPointer(ofocusp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, NULL);
        glVertexAttribPointer(oaefp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*3));
        glVertexAttribPointer(oomegaOmegaincp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*6));
        
        glGenBuffers(1, &data->shader_thick_orbit.thick_orbit_ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data->shader_thick_orbit.thick_orbit_ebo);

        {
            // Multiply by two, because each sample corresponds to a pair of vertices.
            // Mutiply by 3 because there are three indices for one triangle.
            unsigned int num_indices = (data->shader_thick_orbit.anomaly_sample_count_to_use - 1) * 2 * 3;
            unsigned int* indices = malloc(num_indices * sizeof(unsigned int));
            unsigned int currentTriangle = 0;
            for(int vertexIndex = 0; vertexIndex < data->shader_thick_orbit.anomaly_sample_count_to_use * 2; vertexIndex++) {
                unsigned int indicesStart = currentTriangle * 3;

                if(currentTriangle % 2 == 0) {
                    *(indices + indicesStart) = vertexIndex;
                    *(indices + indicesStart + 1) = vertexIndex + 1;
                    *(indices + indicesStart + 2) = vertexIndex + 2;
                    currentTriangle++;
                    continue;
                }

                *(indices + indicesStart) = vertexIndex + 2;
                *(indices + indicesStart + 1) = vertexIndex + 1;
                *(indices + indicesStart + 2) = vertexIndex;
                currentTriangle++;
            }

            // Allocate the generated triangles to the GPU and then deallocate the triangles from the CPU, no longer needed.
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, num_indices * sizeof(unsigned int), indices, GL_STATIC_DRAW);
            free(indices);
        }

        // Set up the trail colors buffer
        glGenBuffers(1, &data->shader_thick_orbit.thick_orbit_cbo);
        glBindBuffer(GL_ARRAY_BUFFER, data->shader_thick_orbit.thick_orbit_cbo);
        
        glEnableVertexAttribArray(trail_colors);
        glVertexAttribPointer(trail_colors, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), NULL);

        // Apply colors
        glBindBuffer(GL_ARRAY_BUFFER, data->shader_thick_orbit.thick_orbit_cbo);
        struct reb_vec3df colors_to_allocate[r->N];
        for(int i = 0; i < r->N; i++) {
            if(i < r->N_orbit_colors) {
                struct reb_vec3df curr_color = r->orbit_colors[i];
                colors_to_allocate[i].x = curr_color.x;
                colors_to_allocate[i].y = curr_color.y;
                colors_to_allocate[i].z = curr_color.z;
                continue;
            }

            // Default to white
            colors_to_allocate[i].x = 1.f;
            colors_to_allocate[i].y = 1.f;
            colors_to_allocate[i].z = 1.f;
        }

        glBufferData(GL_ARRAY_BUFFER, r->N * sizeof(struct reb_vec3df), colors_to_allocate, GL_STATIC_DRAW);

        // Instruct the GPU on what pieces of data attribute to a single rendered instance.
        reb_glVertexAttribDivisor(data->shader_plane.mvp_location, 0); 
        reb_glVertexAttribDivisor(trail_colors, 1);
        reb_glVertexAttribDivisor(ofocusp, 1);
        reb_glVertexAttribDivisor(oaefp, 1);
        reb_glVertexAttribDivisor(oomegaOmegaincp, 1);
        
        // Finished doing allocations and initialization, so unbind the vertex array for clean OpenGL/WebGL state management
        glBindVertexArray(0);
    }

    {
        // PLANE shader
        const char* vertex_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "in vec3 focus;\n"
        "in vec3 aef;\n"
        "in vec3 omegaOmegainc;\n"
        "uniform int vertex_count;\n"
        "uniform mat4 mvp;\n"
        "out float fog;\n"
        "const float M_PI = 3.14159265359;\n"
        "void main() {\n"
        "   float a = aef.x;\n"
        "   float e = aef.y;\n"
        "   float lin = float(gl_VertexID/3)/float(vertex_count/3) + float(gl_VertexID%3)/float(vertex_count/3);\n"
        "   float f = 2.*M_PI*lin;\n"
        "   float theta_max = 0.0;\n"
        "   fog = 1.;\n"
        "   float r;\n"
        "   if (e>1.){\n"
        "       theta_max = acos(-1./e);\n"
        "       lin = 0.5/float(vertex_count/3+1) ;\n"
        "       f = 0.0001-theta_max+1.9998*lin*theta_max;\n"
        "       float rmax = -a*(1.-e*e)/(1. + e*cos(f));\n"
        "       if (gl_VertexID%3==0) { \n"
        "           r = rmax;\n"
        "           f = 0.0; \n"
        "       }else{\n"
        "           lin = float(gl_VertexID/3)/float(vertex_count/3+1) + float(gl_VertexID%3)/float(vertex_count/3+1) - 0.5/float(vertex_count/3+1) ;\n"
        "           f = 0.0001-theta_max+1.9998*lin*theta_max;\n"
        "           r = a*(1.-e*e)/(1. + e*cos(f));\n"
        "       }\n"
        "       fog = 1.-abs(r/rmax);\n"
        "   }else{ \n"
        "       if (gl_VertexID%3==0) { \n"
        "           r = 0.;\n"
        "       }else{\n"
        "           r = a*(1.-e*e)/(1. + e*cos(f));\n"
        "       }\n"
        "   }\n"
        "   float omega = omegaOmegainc.x;\n"
        "   float Omega = omegaOmegainc.y;\n"
        "   float inc = omegaOmegainc.z;\n"
        "   float cO = cos(Omega);\n"
        "   float sO = sin(Omega);\n"
        "   float co = cos(omega);\n"
        "   float so = sin(omega);\n"
        "   float cf = cos(f);\n"
        "   float sf = sin(f);\n"
        "   float ci = cos(inc);\n"
        "   float si = sin(inc);\n"
        "   vec3 pos = vec3(r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci),r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci),+ r*(so*cf+co*sf)*si);\n"
        "   gl_Position = mvp*(vec4(focus+pos, 1.0));\n"
        "   gl_Position.z = 0.;\n" // no clipping
        "}\n";
        const char* fragment_shader =
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;"
        "out vec4 outcolor;\n"
        "in float fog;\n"
        "void main() {\n"
        "  outcolor = vec4(1.,1.,1.,fog*0.3);\n"
        "}\n";

        data->shader_plane.program = loadShader(vertex_shader, fragment_shader);
        data->shader_plane.mvp_location = glGetUniformLocation(data->shader_plane.program, "mvp");
        data->shader_plane.vertex_count_location = glGetUniformLocation(data->shader_plane.program, "vertex_count");

        // Orbit data
        data->shader_plane.vertex_count = 3*200; // higher number = smoother orbits // must be multiple of 3

        // Generate two orbit vao
        glUseProgram(data->shader_plane.program);
        GLuint ofocusp = glGetAttribLocation(data->shader_plane.program,"focus");
        GLuint oaefp = glGetAttribLocation(data->shader_plane.program,"aef");
        GLuint oomegaOmegaincp = glGetAttribLocation(data->shader_plane.program,"omegaOmegainc");

        { // Current
            glGenVertexArrays(1, &data->shader_plane.particle_vao_current);
            glBindVertexArray(data->shader_plane.particle_vao_current);
            glEnableVertexAttribArray(ofocusp);
            glEnableVertexAttribArray(oaefp);
            glEnableVertexAttribArray(oomegaOmegaincp);

            glBindBuffer(GL_ARRAY_BUFFER, data->orbit_buffer_current);
            glVertexAttribPointer(ofocusp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, NULL);
            glVertexAttribPointer(oaefp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*3));
            glVertexAttribPointer(oomegaOmegaincp, 3, GL_FLOAT, GL_FALSE, sizeof(float)*9, (void*)(sizeof(float)*6));

            reb_glVertexAttribDivisor(data->shader_plane.mvp_location, 0); 
            reb_glVertexAttribDivisor(ofocusp, 1);
            reb_glVertexAttribDivisor(oaefp, 1);
            reb_glVertexAttribDivisor(oomegaOmegaincp, 1);
        }

        glBindVertexArray(0);

    }

    { // Set up the post processing objects and buffers
        const char* fullscreen_quad_vertex_shader = 
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "layout (location = 0) in vec3 position;\n"
        "layout (location = 1) in vec2 uvs;\n"
        "out vec2 tex_coords;\n"
        "void main() {\n"
        "   tex_coords = uvs;\n"
        "   gl_Position = vec4(position, 1.f);\n"
        "}";

        // Fragment Shader Source: https://learnopengl.com/Advanced-Lighting/Bloom

        const char* hblur_fragment_shader = 
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;\n"
        "in vec2 tex_coords;\n"
        "out vec4 out_color;\n"

        "uniform sampler2D rendered_texture;\n"

        "const float bloom_size = 2.5f;\n"
        "const int bloom_increments = 5;\n"
        "const float weights[5] = float[] (0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);\n"

        "vec4 apply_horizontal_blur(sampler2D tex, vec2 curr_coord) {\n"
        "   vec4 res = texture(tex, curr_coord) * weights[0];\n"
        "   float tex_offset = 1.f / float(textureSize(tex, 0).x);\n"
        "   for(int i = 1; i < bloom_increments; i++) {\n"
        "       vec2 offset_h = vec2(tex_offset * float(i), 0.f);\n"
        "       res += texture(tex, curr_coord + offset_h) * weights[i];\n"
        "       res += texture(tex, curr_coord - offset_h) * weights[i];\n"
        "   }\n"
        "   return res;\n"
        "}\n"

        "void main() {\n"
        "   out_color = apply_horizontal_blur(rendered_texture, tex_coords);\n"
        "}";

        const char* vblur_and_combine_fragment_shader = 
#ifdef __EMSCRIPTEN__
        "#version 300 es\n"
#else
        "#version 330\n"
#endif
        "precision highp float;\n"
        "in vec2 tex_coords;\n"
        "out vec4 out_color;\n"

        "uniform sampler2D hblur_texture;\n"
        "uniform sampler2D rendered_texture;\n"

        "const float bloom_exposure = 1.125f;\n"
        "const float bloom_size = 2.5f;\n"
        "const int bloom_increments = 5;\n"
        "const float weights[5] = float[] (0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);\n"

        "vec4 apply_vertical_blur(sampler2D tex, vec2 curr_coord) {\n"
        "   vec4 res = texture(tex, curr_coord) * weights[0];\n"
        "   float tex_offset = 1.f / float(textureSize(tex, 0).y);\n"
        "   for(int i = 1; i < bloom_increments; i++) {\n"
        "       float interp = float(i) / float(bloom_increments);\n"
        "       vec2 offset_v = vec2(0.f, tex_offset * float(i));\n"
        "       res += texture(tex, curr_coord + offset_v) * weights[i];\n"
        "       res += texture(tex, curr_coord - offset_v) * weights[i];\n"
        "   }\n"
        "   return res;\n"
        "}\n"

        "void main() {\n"
        "   vec3 unblurred_color = texture(rendered_texture, tex_coords).rgb;\n"
        "   vec3 fullblurred_color = apply_vertical_blur(hblur_texture, tex_coords).rgb;\n"
        "   vec3 added_colors = unblurred_color + fullblurred_color;\n"
        "   vec3 tone_mapped_color = vec3(1.f) - exp(-added_colors * bloom_exposure);\n" // Tone mapping
        "   vec3 bloom_result = pow(tone_mapped_color, vec3(1.f / 1.5f));\n" // Gamma correction with gamma = 1.5
        "   out_color = vec4(bloom_result, 1.f);\n"
        "}";

        // Per line, first three elements are vertex position, last two are UV coordinates
        // (x, y, z, u, v)
        float fullscreen_quad_vertex_data[20] = {
            -1.f, -1.f, 0.f, 0.f, 0.f,
            -1.f, 1.f, 0.f, 0.f, 1.f,
            1.f, 1.f, 0.f, 1.f, 1.f,
            1.f, -1.f, 0.f, 1.f, 0.f
        };

        unsigned int fullscreen_quad_indices[6] = {
            0, 1, 2,
            0, 2, 3
        };

        // The vertex shader is the same for all post processing elements.
        // The only thing that will change between is the fragment shader.
        // This way, all I need to bind is the full screen quad vertex array for rendering through the post processing.
        data->reb_post_processing.horizontal_blur_shader = loadShader(fullscreen_quad_vertex_shader, hblur_fragment_shader);
        data->reb_post_processing.vertical_blur_and_combine_shader = loadShader(fullscreen_quad_vertex_shader, vblur_and_combine_fragment_shader);

        glGenVertexArrays(1, &data->reb_post_processing.fullscreen_quad_vao);
        glBindVertexArray(data->reb_post_processing.fullscreen_quad_vao);
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        glGenBuffers(1, &data->reb_post_processing.fullscreen_quad_vbo);
        glBindBuffer(GL_ARRAY_BUFFER, data->reb_post_processing.fullscreen_quad_vbo);
        glBufferData(GL_ARRAY_BUFFER, 20 * sizeof(float), fullscreen_quad_vertex_data, GL_STATIC_DRAW);

        GLintptr uvs_offset = 3 * sizeof(float);
        GLintptr pos_offset = 0;
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)pos_offset);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)uvs_offset);

        glGenBuffers(1, &data->reb_post_processing.fullscreen_quad_ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, data->reb_post_processing.fullscreen_quad_ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(unsigned int), fullscreen_quad_indices, GL_STATIC_DRAW);
        glBindVertexArray(0);

        // Two frame buffers with a color attachment
        // Frame buffer 1: Have the scene render to this frame buffer which will render to.
        // Frame buffer 2: Will be rendered to with a horizontal blur shader from frame buffer 1
        // Actual window/canvas: Will render with vertical blur and combine with original render from frame buffer 1
        for(int i = 0; i < 2; i++) {
            glGenTextures(1, &data->reb_post_processing.post_processing_color_textures[i]);
            glBindTexture(GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[i]);

            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 700, 700, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
            glBindTexture(GL_TEXTURE_2D, 0);

            // Set up the frame buffer itself
            // If there is a problem with the frame buffer setup, then the keybind to toggle it will be disabled.
            glGenFramebuffers(1, &data->reb_post_processing.post_processing_frame_buffers[i]);
            glBindFramebuffer(GL_FRAMEBUFFER, data->reb_post_processing.post_processing_frame_buffers[i]);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, data->reb_post_processing.post_processing_color_textures[i], 0);
            
            GLenum draw_buffers[1] = { GL_COLOR_ATTACHMENT0 };
            glDrawBuffers(1, draw_buffers);
            
            GLenum frame_buffer_status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
            if(frame_buffer_status != GL_FRAMEBUFFER_COMPLETE) {
                printf("Disabling glow due to frame buffer error. %i", frame_buffer_status);
                data->reb_post_processing.post_processing_enabled = 0;
            } else {
                data->reb_post_processing.post_processing_enabled = 1;
            }

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        // Set up the fragment shader sampler uniforms
        glUseProgram(data->reb_post_processing.fullscreen_quad_shader_program);
        glUniform1i(glGetUniformLocation(data->reb_post_processing.fullscreen_quad_shader_program, "rendered_texture"), 0);

        glUseProgram(data->reb_post_processing.horizontal_blur_shader);
        glUniform1i(glGetUniformLocation(data->reb_post_processing.horizontal_blur_shader, "rendered_texture"), 0);

        glUseProgram(data->reb_post_processing.vertical_blur_and_combine_shader);
        glUniform1i(glGetUniformLocation(data->reb_post_processing.vertical_blur_and_combine_shader, "rendered_texture"), 0);
        glUniform1i(glGetUniformLocation(data->reb_post_processing.vertical_blur_and_combine_shader, "hblur_texture"), 1);

        glUseProgram(0);
    }

    // Main display loop
#ifdef __EMSCRIPTEN__
    // Will return 
    emscripten_request_animation_frame_loop(reb_render_frame_emscripten, r);
#else
    glfwSwapInterval(1);

    while(!glfwWindowShouldClose(window) && r->status<0){
        double t0 = glfwGetTime();
        if (!data->s.pause){
            reb_render_frame(data);
        }

        if (data->take_one_screenshot){
            int cwidth, cheight;
            glfwGetFramebufferSize(data->window, &cwidth, &cheight);
            FILE   *out = fopen("screenshot.tga", "w");
            char   pixel_data[3*cwidth*cheight];
            short  TGAhead[] = {0, 2, 0, 0, 0, 0, cwidth, cheight, 24};

            glReadBuffer(GL_FRONT);
            glReadPixels(0, 0, cwidth, cheight, GL_BGR, GL_UNSIGNED_BYTE, pixel_data);
            fwrite(&TGAhead, sizeof(TGAhead), 1, out);
            fwrite(pixel_data, 3*cwidth*cheight, 1, out);
            fclose(out);
            printf("\nScreenshot saved as 'screenshot.tga',\n");
            data->take_one_screenshot = 0;
        }

        while (glfwGetTime()-t0 < 1.0/120.) { // Maxframerate 120Hz
            usleep(10);
            glfwPollEvents();
        }
    }
    glfwDestroyWindow(window);
    data->window = NULL;
    // Destroy particle buffers
    if (data->N_allocated){
        data->N_allocated = 0;
        data->breadcrumb_N_allocated = 0;
        free(data->particle_data);
        free(data->orbit_data);
        free(data->orbit_colors);
        data->orbit_colors = NULL;
        data->particle_data = NULL;
        data->orbit_data = NULL;
    }
    glfwTerminate();
#endif
}

#endif // OPENGL


