//
//  Solver.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "Particle.h"
#include "Force.h"
#include "BlockSparseMatrix.h"
#include <vector>

struct Vec2fTuple
{
    Vec2f vec1;
    Vec2f vec2;
};

enum class SolverType { Euler, Midpoint, RungeKutta4 };

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                     const std::vector<Force*> & cVector, const SolverType m_SolverType, double* CVector[],
                     double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt);

void EulerSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                 const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                 BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt);

void MidpointSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                    const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                    BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt);

void RungeKutta4thOrderSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                              const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                              BlockSparseMatrix * J, BlockSparseMatrix * JDot,  const double dt);

void ParticleDerivative(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                        const std::vector<Force*> & cVector, std::vector<Vec2fTuple> & dVector);

void Equation11(const std::vector<Particle*> & pVector, const std::vector<Force*> & cVector, double* CVector[],
                double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot);

void ClearForces(const std::vector<Particle*> & pVector);

void CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector);

void ScaleVectorTuples(std::vector<Vec2fTuple> &dVector, const double scaleFactor);