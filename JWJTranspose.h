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
    JWJTranspose(const size_t n, const double W[], BlockSparseMatrix &J);
    ~JWJTranspose();
    void matVecMult(double x[], double r[]) override;
    void matTransVecMult(double x[], double r[]) override;
private:
    const size_t n;
    const double * const W;
    BlockSparseMatrix &J;
    double *m_IntermediateVector;
};

#endif //TUE_SCG_PARTICLE_SYSTEMS_JWJTRANSPOSE_H
