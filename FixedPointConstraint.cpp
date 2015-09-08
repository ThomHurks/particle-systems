#include "FixedPointConstraint.h"
#include <math.h>
#include <assert.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define PI 3.1415926535897932384626433832795

FixedPointConstraint::FixedPointConstraint(
        Particle *p,
        const Vec2f & center,
        BlockSparseMatrix * J,
        BlockSparseMatrix * JDot,
        const int id) :
        // Initialization list:
        Constraint(0, 0), m_p(p), m_center(center),
        m_dC_dp_x(0), m_dC_dp_y(0), m_dCDot_dp_x(0), m_dCDot_dp_y(0)
{
    const int ilength = 1;
    const int jlength = 2;
    // Set up the BSM for the Jacobian:
    double* *dataBlockJ = new double*[ilength * jlength];
    // Todo: assign correct values to the data blocks that are pushed into the BSMs.
    dataBlockJ[0] = &m_dC_dp_x;
    dataBlockJ[1] = &m_dC_dp_y;
    J->AddNewBlock(id, p->m_ID, ilength, jlength, dataBlockJ);
    // Then set up the BSM for the time derivative of the Jacobian:
    double* *dataBlockJDot = new double*[ilength * jlength];
    dataBlockJDot[0] = &m_dCDot_dp_x;
    dataBlockJDot[1] = &m_dCDot_dp_y;
    JDot->AddNewBlock(id, p->m_ID, ilength, jlength, dataBlockJDot);
}

void FixedPointConstraint::draw() const
{
    glBegin(GL_LINE_LOOP);
    glColor3f(1.0,0.5,0.5);
    for (int i=0; i<360; i=i+45)
    {
        double degInRad = i*PI/180;
        glVertex2f(m_center[0] + cos(degInRad) * 0.025, m_center[1] + sin(degInRad) * 0.025);
    }
    glEnd();
}

void FixedPointConstraint::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f relative = m_p->m_Position - m_center;
    m_C = sqrMagnitude(relative) / 2;
    m_CDot = Dot(m_p->m_Velocity, relative);

    Vec2f dC_dp = m_p->m_Position - m_center;
    m_dC_dp_x = dC_dp[0];
    m_dC_dp_y = dC_dp[1];
    m_dCDot_dp_x = m_p->m_Velocity[0];
    m_dCDot_dp_y = m_p->m_Velocity[1];

    assert(!std::isnan(m_C) && std::isfinite(m_C));
    assert(!std::isnan(m_CDot) && std::isfinite(m_CDot));
    assert(!std::isnan(m_dC_dp_x) && std::isfinite(m_dC_dp_x));
    assert(!std::isnan(m_dC_dp_y) && std::isfinite(m_dC_dp_y));
    assert(!std::isnan(m_dCDot_dp_x) && std::isfinite(m_dCDot_dp_x));
    assert(!std::isnan(m_dCDot_dp_y) && std::isfinite(m_dCDot_dp_y));
}