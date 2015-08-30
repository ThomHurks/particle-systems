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
    const int ilength = 1;
    const int jlength = 2;
    // Set up the BSM for the Jacobian:
    double* *dataBlockJ = new double*[ilength * jlength];
    // Todo: assign correct values to the data blocks that are pushed into the BSMs.
    dataBlockJ[0] = &m_C;
    dataBlockJ[1] = &m_C;
    J->AddNewBlock(id, p1->m_ID, ilength, jlength, dataBlockJ);
    J->AddNewBlock(id, p2->m_ID, ilength, jlength, dataBlockJ);
    
    // Then set up the BSM for the time derivative of the Jacobian:
    double* *dataBlockJDot = new double*[ilength * jlength];
    dataBlockJDot[0] = &m_CDot;
    dataBlockJDot[1] = &m_CDot;
    JDot->AddNewBlock(id, p1->m_ID, ilength, jlength, dataBlockJDot);
    JDot->AddNewBlock(id, p2->m_ID, ilength, jlength, dataBlockJDot);
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
