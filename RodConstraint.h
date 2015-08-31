#pragma once

#include "Constraint.h"
#include "BlockSparseMatrix.h"

class RodConstraint : public Constraint
{
public:
    RodConstraint(const Particle *p1, const Particle * p2, const double dist,
                  BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id);
    void draw() const override;
    void ApplyForce(const std::vector<Particle*> & pVector) override;
private:
    const Particle * const m_p1;
    const Particle * const m_p2;
    const double m_distSquared;
};
