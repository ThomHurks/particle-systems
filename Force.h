//
//  Force.h
//  Particle Systems
//
//  Created by Thom Hurks on 21-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#pragma once
#include "Particle.h"
#include <vector>

class Force
{
public:
    Force();
    virtual ~Force() = 0;
    virtual void draw() = 0;
    virtual void ApplyForce(const std::vector<Particle*> & pVector) = 0;
};