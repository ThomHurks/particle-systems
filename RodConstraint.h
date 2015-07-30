#pragma once

#include "Force.h"
#include "BlockSparseMatrix.h"

class RodConstraint : public Force
{
public:
    RodConstraint(const Particle *p1, const Particle * p2, const double dist, double* CVector[],
                  double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id);
    void draw();
    void ApplyForce(const std::vector<Particle*> & pVector);

private:
    const Particle * const m_p1;
    const Particle * const m_p2;
    const double m_distSquared;
    double m_C;
    double m_CDot;
};
