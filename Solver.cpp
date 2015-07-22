#include "Solver.h"

static const float DAMP = 0.98f;

void simulation_step(std::vector<Particle*> pVector, std::vector<Force*> fVector, float dt)
{
    
    
	int ii, size = pVector.size();
	
	for(ii=0; ii<size; ii++)
	{
		pVector[ii]->m_Position += dt * pVector[ii]->m_Velocity;
		pVector[ii]->m_Velocity = DAMP * pVector[ii]->m_Velocity + GetRandomDirection() * 0.005;
	}
}

Vec2f GetRandomDirection()
{
    float randX = (((rand() % 2000) / 1000.f) - 1.f);
    float randY = (((rand() % 2000) / 1000.f) - 1.f);
    return Vec2f(randX, randY);
}

void ParticleDerivative(std::vector<Particle*> pVector, std::vector<Force*> fVector, std::vector<Derivative> dVector)
{
    ClearForces(pVector);
    CalculateForces(pVector, fVector);
    int n = pVector.size();
    int i;
    for(i = 0; i < n; ++i)
    {
        dVector[i]->XDot = pVector[i]->m_Velocity;
        dVector[i]->VDot = pVector[i]->m_AccumulatedForce / pVector[i]->m_Mass;
    }
}

void ClearForces(const std::vector<Particle*> pVector)
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
    int n = pVector.size();
    for(i = 0; i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

