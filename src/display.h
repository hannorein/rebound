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

/**
 * @brief Internal function to check if display update is needed.
 */
void reb_check_for_display_heartbeat(struct reb_simulation* const r);

/**
 * @brief This function initializes OpenGL and starts the run loop.
 * @param data A struct containing all the data needed by the visualization.
 */
void reb_display_init(struct reb_simulation* const r);

void reb_display_init_data(struct reb_simulation* const r);
int reb_display_copy_data(struct reb_simulation* const r);
void reb_display_prepare_data(struct reb_simulation* const r, int orbits);

#endif
