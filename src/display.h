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

struct reb_simulation;

void reb_display_init(struct reb_simulation* const r);

void reb_display_init_data(struct reb_simulation* const r);

#ifdef OPENGL
#ifdef __EMSCRIPTEN__
// Only needed for emscripten client
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
void reb_display_keyboard(GLFWwindow* window, int key, int scancode, int action, int mods);
#endif // __EMSCRIPTEN__
#endif // OPENGL
#endif
