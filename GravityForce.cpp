//
//  GravityForce.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "GravityForce.h"

const Vec2f GravityForce::m_Gravity = Vec2f(0.0, -9.80665);
 
void GravityForce::draw()
{}

void GravityForce::ApplyForce(const std::vector<Particle*> & pVector)
{
    int i;
    int n = pVector.size();
    for(i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce += pVector[i]->m_Mass * m_Gravity;
    }
}