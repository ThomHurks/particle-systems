#include "AngularSpring.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

AngularSpring::AngularSpring(Particle *massPoint, Particle * p1, Particle *p2, double angle, double ks, double kd) :
m_MassPoint(massPoint), m_p1(p1), m_p2(p2), m_angle(angle), m_ks(ks), m_kd(kd)
{
}

void AngularSpring::draw()
{
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 1);
    glVertex2f(m_MassPoint->m_Position[0], m_MassPoint->m_Position[1]);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 1);
    glVertex2f(m_MassPoint->m_Position[0], m_MassPoint->m_Position[1]);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glEnd();
}

void AngularSpring::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f S1 = m_p1->m_Position - m_MassPoint->m_Position;
    Vec2f S2 = m_p2->m_Position - m_MassPoint->m_Position;

    float S1_mag = magnitude(S1);
    float S2_mag = magnitude(S2);

    double currentAngle = acos(Dot(S1, S2) / (S1_mag * S2_mag));
    double angleDelta = (m_angle - currentAngle) / 2;

    Vec2f t1 = RotateAroundPoint(m_MassPoint->m_Position, m_p1->m_Position, -angleDelta);
    Vec2f t2 = RotateAroundPoint(m_MassPoint->m_Position, m_p2->m_Position, angleDelta);

    Vec2f d1 = t1 - m_p1->m_Position;
    Vec2f d2 = t2 - m_p2->m_Position;

    float d1_mag = magnitude(d1);
    float d2_mag = magnitude(d2);

    Vec2f I1_Dot = m_p1->m_Velocity;
    Vec2f I2_Dot = m_p2->m_Velocity;

    Vec2f F1 = -((m_ks * d1_mag) + m_kd * (Dot(I1_Dot, d1) / d1_mag)) * normalized(d1);
    Vec2f F2 = -((m_ks * d2_mag) + m_kd * (Dot(I2_Dot, d2) / d2_mag)) * normalized(d2);

    if (!isnan(F1))
    {  m_p1->m_AccumulatedForce += F1; }

    if (!isnan(F2))
    { m_p2->m_AccumulatedForce += F2; }
}