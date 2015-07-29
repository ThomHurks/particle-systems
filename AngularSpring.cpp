#include "AngularSpring.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdio.h>

AngularSpring::AngularSpring(Particle *p1, Particle * p2, Particle *p3, double angle, double ks, double kd) :
m_p1(p1), m_p2(p2), m_p3(p3), m_angle(angle), m_ks(ks), m_kd(kd)
{
}

void AngularSpring::draw()
{
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 1);
    glVertex2f(m_p1->m_Position[0], m_p1->m_Position[1]);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glEnd();
    
    glBegin(GL_LINES);
    glColor3f(0.3, 0.3, 1);
    glVertex2f(m_p2->m_Position[0], m_p2->m_Position[1]);
    glVertex2f(m_p3->m_Position[0], m_p3->m_Position[1]);
    glEnd();
}

void AngularSpring::ApplyForce(const std::vector<Particle*> & pVector)
{
    Vec2f y1 = m_p2->m_Position-m_p1->m_Position;
    Vec2f y2 = m_p3->m_Position-m_p2->m_Position;
    Vec2f y1dot = m_p2->m_Velocity-m_p1->m_Velocity;
    Vec2f y2dot = m_p3->m_Velocity-m_p2->m_Velocity;
    double y1doty2 = y1*y2;
    double y1mag = magnitude(y1);
    double y2mag = magnitude(y2);
    double y1y2mag = y1mag*y2mag;
    double y1y2mag2 = y1y2mag*y1y2mag;
    double C = y1doty2/y1y2mag-cos(m_angle);
    double Cdotn1=((y1dot*y2)+(y1*y2dot))*y1y2mag;
    double Cdotn2=(y1dot*y1)*(y2*y2)+(y1*y1)*(y2dot*y2);
    double Cdotn3=y1doty2*Cdotn2/y1y2mag;
    double Cdot = (Cdotn1-Cdotn3)/y1y2mag2;
    //printf("  %f,  %f  \n",C,Cdotn1);
    Vec2f temp1 = (1.0/y1y2mag)*y2;
    Vec2f temp2 = (y1doty2/(y1y2mag*y1mag*y1mag))*y1;
    Vec2f deltaCdeltaX1= temp1-temp2;
    
    Vec2f baseF = (-m_ks*C-m_kd*Cdot)*deltaCdeltaX1;
    printf("C:  %f  \n",C);
    printf("Cdot:  %f  \n",Cdot);
    printf("1:  %f,  %f  \n",temp1[0],temp1[1]);
    printf("2:  %f,  %f  \n",temp2[0],temp2[1]);
    
    m_p1->m_AccumulatedForce -= baseF;
    m_p2->m_AccumulatedForce += baseF;
    m_p3->m_AccumulatedForce -= baseF;
}
