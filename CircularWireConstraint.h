#pragma once

#include "Force.h"

class CircularWireConstraint : public Force {
public:
    CircularWireConstraint(Particle *p, const Vec2f & center, const double radius);
    void draw();
    void ApplyForce(const std::vector<Particle*> & pVector);

private:
    static void draw_circle(const Vec2f & vect, float radius);
    Particle * const m_p;
    Vec2f const m_center;
    double const m_radius;
    double const m_radiusSquared;
};
