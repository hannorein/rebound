#ifndef _OUTPUT_H
#define _OUTPUT_H

int output_check(double interval);

void output_ascii(char* filename);
void output_orbit(char* filename);
void output_append_ascii(char* filename);
void output_binary(char* filename);
void output_binary_positions(char* filename);
void output_append_velocity_dispersion(char* filename);
void output_timing();
#ifdef OPENGL
#ifdef LIBPNG
void output_png(char* filename);
#endif
#endif

#endif
