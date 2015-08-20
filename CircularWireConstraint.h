#pragma once

#include "Force.h"
#define PI 3.1415926535897932384626433832795

class CircularWireConstraint : public Force {
public:
    CircularWireConstraint(Particle *p, const Vec2f & center, const double radius);
    void draw() override;
    void ApplyForce(const std::vector<Particle*> & pVector) override;

private:
    static void draw_circle(const Vec2f & vect, const double radius);
    Particle * const m_p;
    Vec2f const m_center;
    double const m_radius;
    double const m_radiusSquared;
};
