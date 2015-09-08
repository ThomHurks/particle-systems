#include "SpringForce.h"
#include <assert.h>
#ifdef __APPLE__
#include <GLUT/glut.h>

#else
#include <GL/glut.h>
#endif

SpringForce::SpringForce(Particle *p1, Particle * p2, double dist, double ks, double kd) :
  m_p1(p1), m_p2(p2), m_dist(dist), m_ks(ks), m_kd(kd) {}

void SpringForce::draw() const
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
    float scale;
    if(l_magnitude == 0)
    {
        scale = m_ks *  m_dist;
    }else
    {
        scale = -(m_ks * (l_magnitude - m_dist) + m_kd * (Dot(lDot, l) / l_magnitude));
    }
    Vec2f f_a = scale * l_unit;
    Vec2f f_b = -f_a;

    assert(!isnan(f_a) && isfinite(f_a));
    assert(!isnan(f_b) && isfinite(f_b));
    
    m_p1->m_AccumulatedForce += f_a;
    m_p2->m_AccumulatedForce += f_b;

    assert(!isnan(m_p1->m_AccumulatedForce) && isfinite(m_p1->m_AccumulatedForce));
    assert(!isnan(m_p2->m_AccumulatedForce) && isfinite(m_p2->m_AccumulatedForce));
}
