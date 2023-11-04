#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "fmemopen.h"
#include "display.h"
#include "input.h"
#include <emscripten/fetch.h>
int first = 1;

void newCallback(struct reb_simulation* r);

void downloadSucceeded(emscripten_fetch_t *fetch) {
    struct reb_simulation* r = fetch->userData;

    FILE* fin = reb_fmemopen((void*)fetch->data, fetch->numBytes, "r");
    enum reb_simulation_binary_error_codes* warnings;
    reb_input_fields(r, fin, warnings);
    fclose(fin);

    if (first){
        r->display_data = calloc(1, sizeof(struct reb_display_data));
        r->display_data->r = r;
        reb_display_init(r); // Will return. Display routines running in animation_loop.
    }
    first = 0;
    emscripten_fetch_close(fetch); // Free data associated with the fetch.

    emscripten_sleep(1000./120.);
    newCallback(r);
}

void downloadFailed(emscripten_fetch_t *fetch) {
    printf("Requesting  %s failed (status code: %d). Server might have shut down.\n", fetch->url, fetch->status);
    emscripten_fetch_close(fetch); // Also free data on failure.
}
void newCallback(struct reb_simulation* r){
    emscripten_fetch_attr_t attr;
    emscripten_fetch_attr_init(&attr);
    attr.userData = r;
    strcpy(attr.requestMethod, "GET");
    attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
    attr.onsuccess = downloadSucceeded;
    attr.onerror = downloadFailed;
    emscripten_fetch(&attr, "http://localhost:1234/simulation");

}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->t = 123.;
    newCallback(r);
    //while(1){
    //sleep(1);
    //}
}
