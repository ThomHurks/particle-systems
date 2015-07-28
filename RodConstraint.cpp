#include "RodConstraint.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

RodConstraint::RodConstraint(Particle *p1, Particle * p2, double dist, BlockSparseMatrix * bsp, int id) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_distSquared(dist * dist), m_C(3) // TEMPORARY TEST, CHANGE TO 0 (TODO)
{
    MatrixBlock block1;
    block1.ci = id;
    block1.pi = p1->m_ID;
    block1.data = &m_C;
    bsp->AddBlock(block1);
    MatrixBlock block2;
    block2.ci = id;
    block2.pi = p2->m_ID;
    block2.data = &m_C;
    bsp->AddBlock(block2);
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
    // Then plug into equation 11 and solve.
}
