//
//  Solver.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#pragma once

#include "Particle.h"
#include "Force.h"
#include "BlockSparseMatrix.h"
#include <vector>
#include "JWJTranspose.h"
#include "linearSolver.h"


class Solver {

    struct Vec2fTuple {
        Vec2f vec1;
        Vec2f vec2;
    };

public:

    enum class SolverType {
        Euler, Midpoint, RungeKutta4
    };

    Solver(const std::vector<Particle *> &pVector, const std::vector<Force *> &fVector,
           const std::vector<Force *> &cVector, double* CVector[],
           double* CDotVector[], BlockSparseMatrix &J, BlockSparseMatrix &JDot);

    void simulation_step(const double dt, SolverType solverType);

    void EulerSolver(const double dt);

    void MidpointSolver(const double dt);

    void RungeKutta4thOrderSolver(const double dt);

    void ParticleDerivative(std::vector<Vec2fTuple> & dVector);

    void Equation11(double ks, double kd);

    void ClearForces(const std::vector<Particle *> &pVector);

    void CalculateForces(const std::vector<Particle *> &pVector, const std::vector<Force *> &fVector);

    void ScaleVectorTuples(std::vector<Vec2fTuple> &dVector, const double scaleFactor);

private:
    const std::vector<Particle *> &m_ParticlesVector;
    const std::vector<Force *> &m_ForcesVector;
    const std::vector<Force *> &m_ConstraintsVector;
    double* *m_CVector;
    double* *m_CDotVector;
    BlockSparseMatrix &m_J;
    BlockSparseMatrix &m_JDot;
};