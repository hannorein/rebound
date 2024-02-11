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
double reconnect_delay = 1.;

void request_frame_from_server(struct reb_simulation* r);

EM_JS(void, reb_hide_console, (int hide), {
    var output = document.getElementById("output");
    if (output){
        if (hide){
        output.style.display = "none";
        }else{
        output.style.display = "block";
        }
    }
    var container = document.getElementById("container");
    if (container){
        if (hide){
            container.style.height = "100%";
        }else{
            container.style.height = "80%";
        }
    }
});

static void request_key_succeeded(emscripten_fetch_t *fetch) {
    emscripten_fetch_close(fetch); // Free data associated with the fetch.
}

void request_frame_failed(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;
    r->display_data->connection_status = -1;
    reconnect_delay *= 1.1;
    printf("Requesting %s failed (status code: %d). Server might have shut down.\n", fetch->url, fetch->status);
    emscripten_fetch_close(fetch); // Also free data on failure.
    
    // Try again after a delay
    emscripten_sleep(1000./120.*reconnect_delay);
    request_frame_from_server(r);
}

void request_key_failed(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;
    r->display_data->connection_status = -1;
    printf("Requesting %s failed (status code: %d). Server might have shut down.\n", fetch->url, fetch->status);
    emscripten_fetch_close(fetch); // Also free data on failure.
}

void send_key(int key){
    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    strcpy(attr.requestMethod, "GET");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = request_key_succeeded;
    attr.onerror = request_key_failed;
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
        send_key(key);
    }
    reb_display_keyboard(window, key, scancode, action, mods);
}


void send_screenshot_succeeded(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;
    emscripten_fetch_close(fetch); // Free data associated with the fetch.
    r->status = REB_STATUS_PAUSED; // Pause until server sends new simulation.
    r->display_data->r_copy->status = REB_STATUS_PAUSED; // Pause until server sends new simulation.
    free(r->display_data->screenshot);
    r->display_data->screenshot = NULL;
    emscripten_sleep(1000./120.);
    request_frame_from_server(r);
}


void send_screenshot_to_server(struct reb_simulation* r){
    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    attr.userData = r;
    strcpy(attr.requestMethod, "POST");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = send_screenshot_succeeded;
    attr.onerror = request_frame_failed;
    attr.requestData = r->display_data->screenshot;
    attr.requestDataSize = strlen(r->display_data->screenshot)+1;

    emscripten_fetch(&attr, "/screenshot");
}


void request_frame_succeeded(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;

    FILE* fin = reb_fmemopen((void*)fetch->data, fetch->numBytes, "r");
    enum reb_simulation_binary_error_codes warnings;
    reb_input_fields(r, fin, &warnings);
    fclose(fin);

    if (first){
        r->display_data = calloc(1, sizeof(struct reb_display_data));
        r->display_data->r = r;
        reb_display_init(r); // Will return. Display routines running in animation_loop.
        glfwSetKeyCallback(r->display_data->window, reb_display_keyboard_passthrough);
        reb_hide_console(1);
    }
    r->display_data->connection_status = 1;
    reconnect_delay = 1.;
    first = 0;
    emscripten_fetch_close(fetch); // Free data associated with the fetch.

    emscripten_sleep(1000./120.);
    request_frame_from_server(r);
}


void request_frame_from_server(struct reb_simulation* r){
    switch (r->status){
        case REB_STATUS_SCREENSHOT:
            // Screenshot not ready yet. Wait. 
            emscripten_sleep(1000./120.);
            request_frame_from_server(r);
            break;
        case REB_STATUS_SCREENSHOT_READY:
            // Screenshot ready
            // Send back
            // resume normal pulls
            send_screenshot_to_server(r);
            break;
        default:
            {
                emscripten_fetch_attr_t attr;
                emscripten_fetch_attr_init(&attr);
                attr.userData = r;
                strcpy(attr.requestMethod, "GET");
                attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
                attr.onsuccess = request_frame_succeeded;
                attr.onerror = request_frame_failed;
                emscripten_fetch(&attr, "/simulation");
            }
            break;
    }


}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->t = 123.;
    request_frame_from_server(r);
    //while(1){
    //sleep(1);
    //}
}
