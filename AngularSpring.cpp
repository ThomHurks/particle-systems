#include "AngularSpring.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

AngularSpring::AngularSpring(Particle *p1, Particle * p2, Particle *p3, double angle, double ks, double kd) :
m_p1(p1), m_p2(p2), m_p3(p3), m_angle(angle), m_ks(ks), m_kd(kd)
{
}

void AngularSpring::draw()
{
    glBegin(GL_LINES);
    glColor3f(0.6, 0.7, 0.8);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glColor3f(0.6, 0.7, 0.8);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glEnd();
}

void AngularSpring::ApplyForce(const std::vector<Particle*> & pVector)
{
    
}
