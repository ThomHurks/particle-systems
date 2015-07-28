#include "MouseSpringForce.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

MouseSpringForce::MouseSpringForce(Particle *p1, double dist, double ks, double kd) :
  m_p(p1), m_dist(dist), m_ks(ks), m_kd(kd) {}

void MouseSpringForce::draw()
{
  glBegin( GL_LINES );
  glColor3f(1.0, 0.7, 0.8);
  glVertex2f( m_loc[0], m_loc[1] );
  glVertex2f( m_p->m_Position[0], m_p->m_Position[1] );
  glEnd();
}

void MouseSpringForce::setMouseLoc(Vec2f loc)
{
    m_loc = loc;
}


void MouseSpringForce::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f l = m_loc - m_p->m_Position;
    float l_magnitude = magnitude(l);
    Vec2f lDot = m_p->m_Velocity;
    Vec2f l_unit = normalized(l);
    
    Vec2f f_a = (m_ks * (l_magnitude - m_dist) + m_kd * ((lDot * l) / l_magnitude)) * l_unit;
    
    m_p->m_AccumulatedForce += f_a;
}

