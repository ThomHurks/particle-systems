#include "Solver.h"

void simulation_step(std::vector<Particle*> pVector, std::vector<Force*> fVector, const float dt)
{
    std::vector<Derivative> dVector(pVector.size()); // Optimize this by reusing.
    ParticleDerivative(pVector, fVector, dVector);
    ScaleDerivativeVector(dVector, dt);
    
	int ii, size = pVector.size();
	
	for(ii=0; ii<size; ii++)
	{
        pVector[ii]->m_Position += dVector[ii].XDot;
        pVector[ii]->m_Velocity += dVector[ii].VDot;
	}
}

void ParticleDerivative(std::vector<Particle*> pVector, std::vector<Force*> fVector, std::vector<Derivative> & dVector)
{
    ClearForces(pVector);
    CalculateForces(pVector, fVector);
    int n = dVector.size();
    int i;
    for(i = 0; i < n; ++i)
    {
        dVector[i].XDot = pVector[i]->m_Velocity;
        dVector[i].VDot = pVector[i]->m_AccumulatedForce / pVector[i]->m_Mass;
    }
}

void ClearForces(std::vector<Particle*> pVector)
{
    int i;
    int n = pVector.size();
    for(i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce = 0;
    }
}

void CalculateForces(std::vector<Particle*> pVector, std::vector<Force*> fVector)
{
    int i;
    int n = fVector.size();
    for(i = 0; i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

void ScaleDerivativeVector(std::vector<Derivative> dVector, const float scaleFactor)
{
    int i;
    int n = dVector.size();
    for(i = 0; i < n; ++i)
    {
        dVector[i].XDot *= scaleFactor;
        dVector[i].VDot *= scaleFactor;
    }
}

