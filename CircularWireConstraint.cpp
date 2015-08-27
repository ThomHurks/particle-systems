#include "CircularWireConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <math.h>

#else
#include <GL/glut.h>
#endif

#define PI 3.1415926535897932384626433832795

CircularWireConstraint::CircularWireConstraint(Particle *p, const Vec2f & center, const double radius) :
        Constraint(0, 0), m_p(p), m_center(center), m_radius(radius), m_radiusSquared(radius * radius) {}

void CircularWireConstraint::draw()
{
    draw_circle(m_center, m_radius);
}

void CircularWireConstraint::ApplyForce(const std::vector<Particle*> & pVector)
{
    m_C = sqrMagnitude(m_p->m_Position - m_center) - m_radiusSquared;
    Vec2f normalAtPoint = normalized(m_p->m_Position - m_center);
    m_CDot = Dot(m_p->m_Velocity, normalAtPoint); // Todo: check if this is correct.
}

void CircularWireConstraint::SetCSlice(double C[])
{
    C[m_p->m_ID] += m_C;
}

void CircularWireConstraint::SetCDotSlice(double CDot[])
{
    CDot[m_p->m_ID] += m_CDot;
}

void CircularWireConstraint::draw_circle(const Vec2f & vect, const double radius)
{
    glBegin(GL_LINE_LOOP);
    glColor3f(0.0,1.0,0.0);
    for (int i=0; i<360; i=i+18)
    {
        double degInRad = i*PI/180;
        glVertex2f(vect[0]+cos(degInRad)*radius,vect[1]+sin(degInRad)*radius);
    }
    glEnd();
}
