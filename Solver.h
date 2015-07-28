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

enum class SolverType { Euler, Midpoint, RungeKutta4 };

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt,const SolverType m_SolverType);

void EulerSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt);

void MidpointSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt);

void RungeKutta4thOrderSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt);

void ParticleDerivative(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, std::vector<Derivative> & dVector);

void ClearForces(const std::vector<Particle*> & pVector);

void CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector);

void ScaleDerivativeVector(std::vector<Derivative> & dVector, const float scaleFactor);