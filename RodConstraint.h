#pragma once

#include "Force.h"
#include "BlockSparseMatrix.h"

class RodConstraint : public Force
{
public:
    RodConstraint(Particle *p1, Particle * p2, double dist, BlockSparseMatrix * bsp, int id);
    void draw();
    void ApplyForce(const std::vector<Particle*> & pVector);

private:
    Particle * const m_p1;
    Particle * const m_p2;
    double const m_dist;
    double const m_distSquared;
    double m_C;
};
