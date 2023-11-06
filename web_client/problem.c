#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <emscripten/fetch.h>
#include "rebound.h"
#include "fmemopen.h"
#include "display.h"
#include "input.h"

int first = 1;

static void request_key_succeeded(emscripten_fetch_t *fetch) {
    emscripten_fetch_close(fetch); // Free data associated with the fetch.
}

void request_failed(emscripten_fetch_t *fetch) {
    printf("Requesting  %s failed (status code: %d). Server might have shut down.\n", fetch->url, fetch->status);
    emscripten_fetch_close(fetch); // Also free data on failure.
}

void send_key(int key){
    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    strcpy(attr.requestMethod, "GET");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = request_key_succeeded;
    attr.onerror = request_failed;
    char buffer[1024];
    sprintf(buffer, "/keyboard/%d", key);
    emscripten_fetch(&attr, buffer);
}

void reb_display_keyboard_passthrough(GLFWwindow* window, int key, int scancode, int action, int mods){
    struct reb_display_data* data = glfwGetWindowUserPointer(window);
    if (!data){
        printf("Error accessing data in reb_display_keyboard\n");
        return;
    }
    if (action==GLFW_PRESS){
        switch(key){
            case 'Q':
            case ' ':
                send_key(key);
                break;
            default:
                break;
        }
    }
    reb_display_keyboard(window, key, scancode, action, mods);
}


void request_frame_from_server(struct reb_simulation* r);

void request_frame_succeeded(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;

    if (first){
    FILE* fin = reb_fmemopen((void*)fetch->data, fetch->numBytes, "r");
    enum reb_simulation_binary_error_codes warnings;
    reb_input_fields(r, fin, &warnings);
    }

    if (first){
        r->display_data = calloc(1, sizeof(struct reb_display_data));
        r->display_data->r = r;
        reb_display_init(r); // Will return. Display routines running in animation_loop.
        glfwSetKeyCallback(r->display_data->window, reb_display_keyboard_passthrough);
    }
    first = 0;
    emscripten_fetch_close(fetch); // Free data associated with the fetch.

    //emscripten_sleep(1000./120.);
    //request_frame_from_server(r);
}


void request_frame_from_server(struct reb_simulation* r){
    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    attr.userData = r;
    strcpy(attr.requestMethod, "GET");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = request_frame_succeeded;
    attr.onerror = request_failed;
    emscripten_fetch(&attr, "/simulation");

}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->t = 123.;
    request_frame_from_server(r);
    //while(1){
    //sleep(1);
    //}
}
