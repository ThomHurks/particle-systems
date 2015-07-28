#include "Solver.h"

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector, const std::vector<Force*> & cVector, const float dt, const SolverType m_SolverType)
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
    int ii, size = pVector.size();
    std::vector<Derivative> dVector1(pVector.size()); // Optimize this by reusing.
    std::vector<Derivative> dVector2(pVector.size()); // Optimize this by reusing.
    std::vector<Derivative> dVector3(pVector.size()); // Optimize this by reusing.
    std::vector<Derivative> dVector4(pVector.size()); // Optimize this by reusing.
    std::vector<Derivative> originals(pVector.size());
    for (ii = 0; ii < size; ii++)
    {
        originals[ii].XDot = pVector[ii]->m_Position;
        originals[ii].VDot = pVector[ii]->m_Velocity;
    }
    //step 1
    ParticleDerivative(pVector, fVector, cVector, dVector1);
    ScaleDerivativeVector(dVector1, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].XDot + dVector1[ii].XDot/2;
        pVector[ii]->m_Velocity = originals[ii].VDot + dVector1[ii].VDot/2;
    }
    //step 2
    ParticleDerivative(pVector, fVector, cVector, dVector2);
    ScaleDerivativeVector(dVector2, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].XDot + dVector2[ii].XDot/2;
        pVector[ii]->m_Velocity = originals[ii].VDot + dVector2[ii].VDot/2;
    }
    //step 3
    ParticleDerivative(pVector, fVector, cVector, dVector3);
    ScaleDerivativeVector(dVector3, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].XDot + dVector3[ii].XDot;
        pVector[ii]->m_Velocity = originals[ii].VDot + dVector3[ii].VDot;
    }
    //step 4
    ParticleDerivative(pVector, fVector, cVector, dVector4);
    ScaleDerivativeVector(dVector4, dt);
    //final positions & velocities
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].XDot + (dVector1[ii].XDot+2*dVector2[ii].XDot+2*dVector3[ii].XDot+dVector4[ii].XDot)/6;
        pVector[ii]->m_Velocity = originals[ii].VDot + (dVector1[ii].VDot+2*dVector2[ii].VDot+2*dVector3[ii].VDot+dVector4[ii].VDot)/6;
    }
    
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

