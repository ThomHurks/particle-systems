//
//  BlockSparseMatrix.h
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#pragma once
#include "linearSolver.h"
#include <vector>

struct MatrixBlock {
    const int ci; // constraint index (i)
    const int pi; // particle index(j)
    const int ilength;//width of block (always 1 in our case)
    const int jlength;//height of block (always 2 in our case)

    const double * const data;

    MatrixBlock(const int ci, const int pi, const int ilength, const int jlength, const double * const data) :
    ci(ci), pi(pi), ilength(ilength), jlength(jlength), data(data) {}
};

class BlockSparseMatrix : public implicitMatrixWithTrans
{
public:
    void matVecMult(double x[], double r[]) override;
    void matTransVecMult(double x[], double r[]) override;
    void AddNewBlock(const int ci, const int pi, const int ilength, const int jlength, const double * const data);
private:
    std::vector<MatrixBlock> m_Matrix;
};