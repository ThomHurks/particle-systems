#pragma once

#include "Constraint.h"
#include "BlockSparseMatrix.h"

#define PI 3.1415926535897932384626433832795

class FixedPointConstraint : public Constraint {
public:
    FixedPointConstraint(Particle *p, const Vec2f & center,
                           BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id);
    void draw() const override;
    void ApplyForce(const std::vector<Particle*> & pVector) override;
private:
    Particle * const m_p;
    Vec2f const m_center;
    double m_dC_dp_x;
    double m_dC_dp_y;
    double m_dCDot_dp_x;
    double m_dCDot_dp_y;
};
