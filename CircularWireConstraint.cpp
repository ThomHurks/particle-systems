#include "CircularWireConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <math.h>

#else
#include <GL/glut.h>
#endif

#define PI 3.1415926535897932384626433832795

CircularWireConstraint::CircularWireConstraint(
Particle *p,
        const Vec2f & center,
        const double radius,
        BlockSparseMatrix * J,
        BlockSparseMatrix * JDot,
        const int id) :
        Constraint(0, 0), m_p(p), m_center(center), m_radius(radius), m_radiusSquared(radius * radius)
{
    const int ilength = 1;
    const int jlength = 2;
    // Set up the BSM for the Jacobian:
    double* *dataBlockJ = new double*[ilength * jlength];
    // Todo: assign correct values to the data blocks that are pushed into the BSMs.
    dataBlockJ[0] = &m_C;
    dataBlockJ[1] = &m_C;
    J->AddNewBlock(id, p->m_ID, ilength, jlength, dataBlockJ);
    // Then set up the BSM for the time derivative of the Jacobian:
    double* *dataBlockJDot = new double*[ilength * jlength];
    dataBlockJDot[0] = &m_CDot;
    dataBlockJDot[1] = &m_CDot;
    JDot->AddNewBlock(id, p->m_ID, ilength, jlength, dataBlockJDot);
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
    Vec2f relative = m_p->m_Position - m_center;
    m_C = (sqrMagnitude(relative) - m_radiusSquared) / 2;
    m_CDot = Dot(m_p->m_Velocity, relative);
    //m_dCdx = relative;
    //m_dCdotdx = m_p->m_Velocity;
}

void CircularWireConstraint::SetCSlice(double C[]) const
{
    C[m_p->m_ID] += m_C;
}

void CircularWireConstraint::SetCDotSlice(double CDot[]) const
{
    CDot[m_p->m_ID] += m_CDot;
}