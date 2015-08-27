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

// Abstract class
class Force
{
public:
    virtual ~Force() {};
    virtual void draw() const = 0;
    virtual void ApplyForce(const std::vector<Particle*> & pVector) = 0;
};