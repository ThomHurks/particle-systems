#pragma once

#include "Particle.h"
#include "Force.h"

class AngularSpring : public Force {
    
 public:
  AngularSpring(Particle *p1, Particle * p2,Particle * p3, double angle, double ks, double kd);
  void draw();
  void ApplyForce(const std::vector<Particle*> & pVector);

 private:
  Particle * const m_p1;   // particle 1
  Particle * const m_p2;   // particle 2
  Particle * const m_p3;   // particle 3 
  double const m_angle;     // rest angle
  double const m_ks, m_kd; // spring strength constants
};

