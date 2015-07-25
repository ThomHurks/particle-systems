#include "Solver.h"

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt)
{
    SolverType m_SolverType = SolverType::Euler; // Use this to set the solver type.
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

void EulerSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt)
{
    int ii, size = pVector.size();
    std::vector<Derivative> dVector(pVector.size()); // Optimize this by reusing.
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleDerivativeVector(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position += dVector[ii].XDot;
        pVector[ii]->m_Velocity += dVector[ii].VDot;
    }
}

void MidpointSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt)
{
    int ii, size = pVector.size();
    std::vector<Derivative> dVector(pVector.size()); // Optimize this by reusing.
    std::vector<Derivative> originals(pVector.size());
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleDerivativeVector(dVector, dt / 2); //notice the /2
    for (ii = 0; ii < size; ii++)
    {   //Make half an euler step, save original positions and velocities.
        originals[ii].XDot = pVector[ii]->m_Position;
        originals[ii].VDot = pVector[ii]->m_Velocity;
        pVector[ii]->m_Position += dVector[ii].XDot;
        pVector[ii]->m_Velocity += dVector[ii].VDot;
    }
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleDerivativeVector(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].XDot + dVector[ii].XDot;
        pVector[ii]->m_Velocity = originals[ii].VDot + dVector[ii].VDot;
    }
}

void RungeKutta4thOrderSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt)
{
    
}

void ParticleDerivative(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, std::vector<Derivative> & dVector)
{
    ClearForces(pVector);
    // First calculate forces:
    CalculateForces(pVector, fVector);
    // Then calculate constraint forces:
    CalculateForces(pVector, cVector);
    // Then derivatives:
    int n = dVector.size();
    int i;
    for (i = 0; i < n; ++i)
    {
        dVector[i].XDot = pVector[i]->m_Velocity;
        dVector[i].VDot = pVector[i]->m_AccumulatedForce / pVector[i]->m_Mass;
    }
}

void ClearForces(const std::vector<Particle*> & pVector)
{
    int i;
    int n = pVector.size();
    for (i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce = 0;
    }
}

void CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector)
{
    int i;
    int n = fVector.size();
    for (i = 0; i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

void ScaleDerivativeVector(std::vector<Derivative> & dVector, const float scaleFactor)
{
    int i;
    int n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].XDot *= scaleFactor;
        dVector[i].VDot *= scaleFactor;
    }
}

