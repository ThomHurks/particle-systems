//
//  GravityForce.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "Force.h"

class GravityForce : public Force
{
public:
    void draw() const override;
    void ApplyForce(const std::vector<Particle*> & pVector) override;
private:
    static const Vec2f m_Gravity;
};