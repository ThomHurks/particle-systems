#include "RodConstraint.h"
#include <assert.h>
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
    // Todo: assign correct values to the data blocks that are pushed into the BSMs.
    double* *dataBlockJp1 = new double*[ilength * jlength];
    dataBlockJp1[0] = &m_dC_dp1_x;//should be x of m_p1.position-m_p2.position
    dataBlockJp1[1] = &m_dC_dp1_y;//should be y of m_p1.position-m_p2.position
    
    double* *dataBlockJp2 = new double*[ilength * jlength];
    dataBlockJp2[0] = &m_dC_dp2_x;//should be x of m_p2.position-m_p1.position
    dataBlockJp2[1] = &m_dC_dp2_y;//should be y of m_p2.position-m_p1.position
    
    J->AddNewBlock(id, p1->m_ID, ilength, jlength, dataBlockJp1);
    J->AddNewBlock(id, p2->m_ID, ilength, jlength, dataBlockJp2);
    
    // Then set up the BSM for the time derivative of the Jacobian:
    double* *dataBlockJDotp1 = new double*[ilength * jlength];
    dataBlockJDotp1[0] = &m_dCDot_dp1_x;//should be x of m_p1.velocity-m_p2.velocity
    dataBlockJDotp1[1] = &m_dCDot_dp1_y;//should be y of m_p1.velocity-m_p2.velocity
    
    double* *dataBlockJDotp2 = new double*[ilength * jlength];
    dataBlockJDotp2[0] = &m_dCDot_dp2_x;//should be x of m_p2.velocity-m_p1.velocity
    dataBlockJDotp2[1] = &m_dCDot_dp2_y;//should be y of m_p2.velocity-m_p1.velocity
    
    JDot->AddNewBlock(id, p1->m_ID, ilength, jlength, dataBlockJDotp1);
    JDot->AddNewBlock(id, p2->m_ID, ilength, jlength, dataBlockJDotp2);
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
    Vec2f dp = m_p1->m_Position - m_p2->m_Position;
    Vec2f dv = m_p1->m_Velocity - m_p2->m_Velocity;
    m_C = (sqrMagnitude(dp) - m_distSquared)/2;
    m_CDot = Dot(dv,dp);
    
    m_dC_dp1_x=dp[0];//should be x of m_p1.position-m_p2.position
    m_dC_dp1_y=dp[1];//should be y of m_p1.position-m_p2.position
    
    m_dC_dp2_x=-dp[0];//should be x of m_p2.position-m_p1.position
    m_dC_dp2_y=-dp[1];//should be y of m_p2.position-m_p1.position
    
    m_dCDot_dp1_x=dv[0];//should be x of m_p1.velocity-m_p2.velocity
    m_dCDot_dp1_y=dv[1];//should be y of m_p1.velocity-m_p2.velocity
    
    m_dCDot_dp2_x=-dv[0];//should be x of m_p2.velocity-m_p1.velocity
    m_dCDot_dp2_y=-dv[1];//should be y of m_p2.velocity-m_p1.velocity

    assert(!std::isnan(m_C) && std::isfinite(m_C));
    assert(!std::isnan(m_CDot) && std::isfinite(m_CDot));
    assert(!std::isnan(m_dC_dp1_x) && std::isfinite(m_dC_dp1_x));
    assert(!std::isnan(m_dC_dp1_y) && std::isfinite(m_dC_dp1_y));
    assert(!std::isnan(m_dC_dp2_x) && std::isfinite(m_dC_dp2_x));
    assert(!std::isnan(m_dC_dp2_y) && std::isfinite(m_dC_dp2_y));
    assert(!std::isnan(m_dCDot_dp1_x) && std::isfinite(m_dCDot_dp1_x));
    assert(!std::isnan(m_dCDot_dp1_y) && std::isfinite(m_dCDot_dp1_y));
    assert(!std::isnan(m_dCDot_dp2_x) && std::isfinite(m_dCDot_dp2_x));
    assert(!std::isnan(m_dCDot_dp2_y) && std::isfinite(m_dCDot_dp2_y));
}
