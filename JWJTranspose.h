//
// Created by Thom Hurks on 20-08-15.
//

#ifndef TUE_SCG_PARTICLE_SYSTEMS_JWJTRANSPOSE_H
#define TUE_SCG_PARTICLE_SYSTEMS_JWJTRANSPOSE_H

#pragma once

#include "linearSolver.h"
#include "BlockSparseMatrix.h"

class JWJTranspose : public implicitMatrixWithTrans
{
public:
    JWJTranspose(double W[], BlockSparseMatrix *J);
    void matVecMult(double x[], double r[]) override;
    void matTransVecMult(double x[], double r[]) override;
private:
    double *W;
    BlockSparseMatrix *J;
};


#endif //TUE_SCG_PARTICLE_SYSTEMS_JWJTRANSPOSE_H
