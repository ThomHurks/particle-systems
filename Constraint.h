//
// Created by Thom Hurks on 24-08-15.
//

#ifndef TUE_SCG_PARTICLE_SYSTEMS_CONSTRAINT_H
#define TUE_SCG_PARTICLE_SYSTEMS_CONSTRAINT_H


#include "Force.h"

// Abstract class
class Constraint : public Force {
public:
    Constraint(const double C, const double CDot) : m_C(C), m_CDot(CDot) { }
    virtual ~Constraint() {};
    // Getters:
    double GetC() const { return m_C; }
    double GetCDot() const { return m_CDot; }
protected:
    double m_C;
    double m_CDot;
};


#endif //TUE_SCG_PARTICLE_SYSTEMS_CONSTRAINT_H
