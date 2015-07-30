#include "Solver.h"

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const double dt, const SolverType m_SolverType)
{
    switch (m_SolverType)
    {
        case SolverType::Euler:
            EulerSolver(pVector, fVector, cVector, dt);
            break;
        case SolverType::Midpoint:
            MidpointSolver(pVector, fVector, cVector, dt);
            break;
        case SolverType::RungeKutta4:
            RungeKutta4thOrderSolver(pVector, fVector, cVector, dt);
            break;
    }
}

void EulerSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position += dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
}

void MidpointSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt / 2); //notice the /2
    for (ii = 0; ii < size; ii++)
    {   //Make half an euler step, save original positions and velocities.
        originals[ii].vec1 = pVector[ii]->m_Position;
        originals[ii].vec2 = pVector[ii]->m_Velocity;
        pVector[ii]->m_Position += dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector[ii].vec2; // VDot
    }
}

void RungeKutta4thOrderSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                              const std::vector<Force*> & cVector, const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector1(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector2(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector3(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector4(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    for (ii = 0; ii < size; ii++)
    {
        originals[ii].vec1 = pVector[ii]->m_Position;
        originals[ii].vec2 = pVector[ii]->m_Velocity;
    }
    //step 1
    ParticleDerivative(pVector, fVector, cVector, dVector1);
    ScaleVectorTuples(dVector1, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector1[ii].vec1 /2; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector1[ii].vec2 /2; // VDot
    }
    //step 2
    ParticleDerivative(pVector, fVector, cVector, dVector2);
    ScaleVectorTuples(dVector2, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector2[ii].vec1 /2; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector2[ii].vec2 /2; // VDot
    }
    //step 3
    ParticleDerivative(pVector, fVector, cVector, dVector3);
    ScaleVectorTuples(dVector3, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector3[ii].vec1; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector3[ii].vec2; // VDot
    }
    //step 4
    ParticleDerivative(pVector, fVector, cVector, dVector4);
    ScaleVectorTuples(dVector4, dt);
    //final positions & velocities
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 +
                (dVector1[ii].vec1 +2*dVector2[ii].vec1 +2*dVector3[ii].vec1 +dVector4[ii].vec1)/6; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 +
                (dVector1[ii].vec2 +2*dVector2[ii].vec2 +2*dVector3[ii].vec2 +dVector4[ii].vec2)/6; // VDot
    }
    
}

void ParticleDerivative(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, std::vector<Vec2fTuple> & dVector)
{
    ClearForces(pVector);
    // First calculate forces:
    CalculateForces(pVector, fVector);
    // Then calculate constraint forces:
    CalculateForces(pVector, cVector);
    // Then derivatives:
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 = pVector[i]->m_Velocity; // XDot
        dVector[i].vec2 = pVector[i]->m_AccumulatedForce / pVector[i]->m_Mass; // VDot
    }
}

// Attempt at implementing equation 11, so far only set up some variables:
void Equation11(const std::vector<Particle*> & pVector, const std::vector<Force*> & cVector)
{
    size_t i, n = pVector.size();
    Vec2f * q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { q[i] = pVector[i]->m_Position; }
    double * M = new double[n];
    for (i = 0; i < n; ++i)
    { M[i] = pVector[i]->m_Mass; }
    Vec2f * Q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { Q[i] = pVector[i]->m_AccumulatedForce; }
    double * W = new double[n];
    for (i = 0; i < n; ++i)
    { W[i] = 1 / pVector[i]->m_Mass; }
}

void ClearForces(const std::vector<Particle*> & pVector)
{
    size_t i, n = pVector.size();
    for (i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce = 0;
    }
}

void CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector)
{
    size_t i, n = fVector.size();
    for (i = 0; i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

void ScaleVectorTuples(std::vector<Vec2fTuple> &dVector, const double scaleFactor)
{
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 *= scaleFactor;
        dVector[i].vec2 *= scaleFactor;
    }
}

