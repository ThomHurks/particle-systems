#pragma once

#include <gfx/vec2.h>

class Particle
{
public:

    Particle(const Vec2f & ConstructPos, int id);
	virtual ~Particle(void);

	void reset();
	void draw() const;

	const Vec2f m_ConstructPos;
	Vec2f m_Position;
	Vec2f m_Velocity;
    Vec2f m_AccumulatedForce;
    float m_Mass;
	int m_ID;
};
