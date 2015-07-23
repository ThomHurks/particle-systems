#include "SpringForce.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

SpringForce::SpringForce(Particle *p1, Particle * p2, double dist, double ks, double kd) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {}

void SpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p1->m_Position[0], m_p1->m_Position[1] );
  glColor3f(0.6, 0.7, 0.8);
  glVertex2f( m_p2->m_Position[0], m_p2->m_Position[1] );
  glEnd();
}

void SpringForce::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f l = m_p1->m_Position - m_p2->m_Position;
    float l_magnitude = magnitude(l);
    Vec2f lDot = m_p1->m_Velocity - m_p2->m_Velocity;
    Vec2f l_unit = normalized(l);
    
    Vec2f f_a = -(m_ks * (l_magnitude - m_dist) + m_kd * ((lDot * l) / l_magnitude)) * l_unit;
    Vec2f f_b = -f_a;
    
    m_p1->m_AccumulatedForce += f_a;
    m_p2->m_AccumulatedForce += f_b;
}
