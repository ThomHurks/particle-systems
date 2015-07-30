#include "RodConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

RodConstraint::RodConstraint(const Particle *p1, const Particle * p2, const double dist, double* CVector[],
                             double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot, const int id) :
  m_p1(p1), m_p2(p2), m_distSquared(dist * dist), m_C(0), m_CDot(0)
{
    // Todo: improve the code duplication in this constructor.
    // Setup C and CDot
    CVector[id] = &m_C;
    CDotVector[id] = &m_CDot;
    // First set up the BSM for the Jacobian:
    MatrixBlock block1;
    block1.ci = id;
    block1.pi = p1->m_ID;
    block1.data = &m_C;
    J->AddBlock(block1);
    MatrixBlock block2;
    block2.ci = id;
    block2.pi = p2->m_ID;
    block2.data = &m_C;
    J->AddBlock(block2);
    // Then set up the BSM for the time derivative of the Jacobian:
    MatrixBlock block3;
    block3.ci = id;
    block3.pi = p1->m_ID;
    block3.data = &m_CDot;
    JDot->AddBlock(block3);
    MatrixBlock block4;
    block4.ci = id;
    block4.pi = p2->m_ID;
    block4.data = &m_CDot;
    JDot->AddBlock(block4);
}

void RodConstraint::draw()
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
    m_CDot = 0; // Todo: Implement this.
    // Then plug into equation 11 and solve (TODO)
}
