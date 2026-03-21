/**
 * @file 	particle.h
 * @brief 	reb_particle structure and main particle routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#ifndef _PARTICLE_H
#define _PARTICLE_H
struct reb_simulation;
struct reb_particle;

// Hash table for names.
// Open Hashing, linked list.
struct reb_name_hash_item {
    size_t index;
    struct reb_name_hash_item* next;
};

// Returns the index of the rootbox for the particle based on its position.
int reb_get_rootbox_for_particle(const struct reb_simulation* const r, struct reb_particle pt);

// Returns 1 if a testparticle of type 0 has a finite mass (it should have mass=0).
int reb_particle_check_testparticles(struct reb_simulation* const r);
#endif // _PARTICLE_H
