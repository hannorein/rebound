#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "fmemopen.h"
#include "display.h"
#include "input.h"
#include <emscripten/fetch.h>


void downloadSucceeded(emscripten_fetch_t *fetch) {
    printf("Finished downloading %llu bytes from URL %s.\n", fetch->numBytes, fetch->url);
    struct reb_simulation* r = fetch->userData;
    // The data is now available at fetch->data[0] through fetch->data[fetch->numBytes-1];

    FILE* fin = reb_fmemopen(fetch->data, fetch->numBytes, "r");
    enum reb_simulation_binary_error_codes* warnings;
    reb_input_fields(r, fin, warnings);
    fclose(fin);

    printf("t=%f\n",r->t);

    emscripten_fetch_close(fetch); // Free data associated with the fetch.
}

void downloadFailed(emscripten_fetch_t *fetch) {
    printf("Downloading %s failed, HTTP failure status code: %d.\n", fetch->url, fetch->status);
    emscripten_fetch_close(fetch); // Also free data on failure.
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->t = 123.;
                    reb_display_init_data(r);
                    r->display_data->opengl_enabled = 1;
                    reb_display_init(r); // Will return. Display routines running in animation_loop.
    //while(1){
      emscripten_fetch_attr_t attr;
      emscripten_fetch_attr_init(&attr);
      attr.userData = r;
      strcpy(attr.requestMethod, "GET");
      attr.attributes = EMSCRIPTEN_FETCH_LOAD_TO_MEMORY;
      attr.onsuccess = downloadSucceeded;
      attr.onerror = downloadFailed;
      printf("fetch now\n");
      emscripten_fetch(&attr, "http://localhost:1234/simulation");
      printf("Sleep.\n");
      //sleep(1);
  //}
}
