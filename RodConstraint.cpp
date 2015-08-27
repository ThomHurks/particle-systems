#include "RodConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

RodConstraint::RodConstraint(const Particle *p1, const Particle * p2, const double dist,
                                          BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id) :
  Constraint(0, 0), m_p1(p1), m_p2(p2), m_distSquared(dist * dist)
{
    // Set up the BSM for the Jacobian:
    J->AddNewBlock(id, p1->m_ID, &m_C);
    J->AddNewBlock(id, p2->m_ID, &m_C);
    // Then set up the BSM for the time derivative of the Jacobian:
    JDot->AddNewBlock(id, p1->m_ID, &m_CDot);
    JDot->AddNewBlock(id, p2->m_ID, &m_CDot);
}

void RodConstraint::draw() const
{
  glBegin( GL_LINES );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.8, 0.7, 0.6);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}

void RodConstraint::ApplyForce(const std::vector<Particle*> & pVector)
{
    m_C = sqrMagnitude(m_p1->m_Position - m_p2->m_Position) - m_distSquared;
    m_CDot = sqrMagnitude(m_p1->m_Velocity - m_p2->m_Velocity); // Todo: check if this is correct.
}

void RodConstraint::SetCSlice(double C[]) const
{
    C[m_p1->m_ID] += m_C;
    C[m_p2->m_ID] += m_C;
}

void RodConstraint::SetCDotSlice(double CDot[]) const
{
    CDot[m_p1->m_ID] += m_CDot;
    CDot[m_p2->m_ID] += m_CDot;
}
