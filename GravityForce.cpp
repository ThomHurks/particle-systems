//
//  GravityForce.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "GravityForce.h"

const Vec2f GravityForce::m_Gravity = Vec2f(0.0, -6.67384);

void GravityForce::draw()
{}

void GravityForce::ApplyForce(std::vector<Particle*> pVector)
{
    int i;
    int n = pVector.size();
    for(i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce += m_Gravity;
    }
}