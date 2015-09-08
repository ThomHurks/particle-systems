#include "MouseSpringForce.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <assert.h>

#else
#include <GL/glut.h>
#endif

MouseSpringForce::MouseSpringForce(Particle *p, const double dist, const double ks, const double kd,
                                   const float mousePosition[2]) :
  m_p(p), m_dist(dist), m_ks(ks), m_kd(kd), m_loc(mousePosition) {}

void MouseSpringForce::draw() const
{
  glBegin( GL_LINES );
  glColor3f(1.0, 0.7, 0.8);
  glVertex2f( m_loc[0], m_loc[1] );
  glVertex2f( m_p->m_Position[0], m_p->m_Position[1] );
  glEnd();
}

void MouseSpringForce::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f l = Vec2f(m_loc) - m_p->m_Position;
    float l_magnitude = magnitude(l);
    Vec2f lDot = m_p->m_Velocity;
    Vec2f l_unit = normalized(l);
    
    Vec2f f_a = (m_ks * (l_magnitude - m_dist) + m_kd * (Dot(lDot, l) / l_magnitude)) * l_unit;

    assert(!isnan(f_a) && isfinite(f_a));
    m_p->m_AccumulatedForce += f_a;
    assert(!isnan(m_p->m_AccumulatedForce) && isfinite(m_p->m_AccumulatedForce));
}