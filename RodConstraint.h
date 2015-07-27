#pragma once

#include "Force.h"

class RodConstraint : public Force
{
public:
    RodConstraint(Particle *p1, Particle * p2, double dist);
    void draw();
    void ApplyForce(const std::vector<Particle*> & pVector);

private:
    Particle * const m_p1;
    Particle * const m_p2;
    double const m_dist;
    double const m_distSquared;
};
