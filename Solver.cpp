#include "Solver.h"

static const float DAMP = 0.98f;

void simulation_step(std::vector<Particle*> pVector, float dt)
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

