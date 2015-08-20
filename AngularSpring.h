#pragma once

#include "Particle.h"
#include "Force.h"

class AngularSpring : public Force {
    
 public:
  AngularSpring(Particle *p1, Particle * p2,Particle * p3, double angle, double ks, double kd);
  void draw() override;
  void ApplyForce(const std::vector<Particle*> & pVector) override;

 private:
  Particle * const m_MassPoint; // The mass point
  Particle * const m_p1;        // particle 1
  Particle * const m_p2;        // particle 2
  double const m_angle;         // cosine of rest angle
  double const m_ks, m_kd;      // spring strength constants
};