/**
 * @file 	rotations.c
 * @brief 	Tools for rotations and quaternions.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2022 Hanno Rein
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
#include <stdlib.h>
#include <math.h>
#include "rebound.h"


struct reb_quat reb_quat_identity(){
    struct reb_quat q = {.ix = 0.0, .iy = 0.0, .iz = 0.0, .r = 1.0 };
    return q;
}
