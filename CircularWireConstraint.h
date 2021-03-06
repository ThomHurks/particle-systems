#pragma once

#include "Constraint.h"
#include "BlockSparseMatrix.h"

#define PI 3.1415926535897932384626433832795

class CircularWireConstraint : public Constraint {
public:
    CircularWireConstraint(Particle *p, const Vec2f & center, const double radius,
                           BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id);
    void draw() const override;
    void ApplyForce(const std::vector<Particle*> & pVector) override;
private:
    Particle * const m_p;
    Vec2f const m_center;
    double const m_radius;
    double const m_radiusSquared;
    double m_dC_dp_x;
    double m_dC_dp_y;
    double m_dCDot_dp_x;
    double m_dCDot_dp_y;
};
