#include "CircularWireConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <math.h>

#else
#include <GL/glut.h>
#endif

#define PI 3.1415926535897932384626433832795

CircularWireConstraint::CircularWireConstraint(Particle *p, const Vec2f & center, const double radius,
                                               BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id) :
        Constraint(0, 0), m_p(p), m_center(center), m_radius(radius), m_radiusSquared(radius * radius)
{
    // Set up the BSM for the Jacobian:
    J->AddNewBlock(id, p->m_ID, &m_C);
    // Then set up the BSM for the time derivative of the Jacobian:
    JDot->AddNewBlock(id, p->m_ID, &m_CDot);
}

void CircularWireConstraint::draw() const
{
    glBegin(GL_LINE_LOOP);
    glColor3f(0.0,1.0,0.0);
    for (int i=0; i<360; i=i+18)
    {
        double degInRad = i*PI/180;
        glVertex2f(m_center[0] + cos(degInRad) * m_radius, m_center[1] + sin(degInRad) * m_radius);
    }
    glEnd();
}

void CircularWireConstraint::ApplyForce(const std::vector<Particle*> & pVector)
{
    m_C = sqrMagnitude(m_p->m_Position - m_center) - m_radiusSquared;
    Vec2f normalAtPoint = normalized(m_p->m_Position - m_center);
    m_CDot = Dot(m_p->m_Velocity, normalAtPoint); // Todo: check if this is correct.
}

void CircularWireConstraint::SetCSlice(double C[]) const
{
    C[m_p->m_ID] += m_C;
}

void CircularWireConstraint::SetCDotSlice(double CDot[]) const
{
    CDot[m_p->m_ID] += m_CDot;
}