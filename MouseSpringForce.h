/* 
 * File:   MouseSpringForce.h
 * Author: Helmond
 *
 * Created on 28 juli 2015, 16:42
 */

#ifndef MOUSESPRINGFORCE_H
#define	MOUSESPRINGFORCE_H
#endif	/* MOUSESPRINGFORCE_H */

#pragma once

#include "Particle.h"
#include "Force.h"

class MouseSpringForce : public Force {
    
 public:
  MouseSpringForce(Particle *p, const double dist, const double ks, const double kd, const float mousePosition[2]);
  void draw() override;
  void ApplyForce(const std::vector<Particle*> & pVector) override;

 private:
  Particle * const m_p;
  double const m_dist;      // rest length
  double const m_ks, m_kd;  // spring strength constants
  const float * const m_loc;// current mouse location
};