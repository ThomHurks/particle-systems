//
//  Solver.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "Particle.h"
#include <vector>

void simulation_step(std::vector<Particle*> pVector, float dt);

Vec2f GetRandomDirection();