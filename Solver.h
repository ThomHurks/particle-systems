//
//  Solver.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "Particle.h"
#include "Force.h"
#include <vector>

struct Derivative {
    Vec2f XDot;
    Vec2f VDot;
};

void simulation_step(std::vector<Particle*> pVector, std::vector<Force*> fVector, const float dt);

void ParticleDerivative(std::vector<Particle*> pVector, std::vector<Force*> fVector, std::vector<Derivative> & dVector);

void ClearForces(std::vector<Particle*> pVector);

void CalculateForces(std::vector<Particle*> pVector, std::vector<Force*> fVector);

void ScaleDerivativeVector(std::vector<Derivative> dVector, const float scaleFactor);