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
    const int ci; // constraint index
    const int pi; // particle index

    double * const data;

    MatrixBlock(const int ci, const int pi, double * const data) : ci(ci), pi(pi), data(data) {}
};

class BlockSparseMatrix : public implicitMatrixWithTrans
{
public:
    void matVecMult(double x[], double r[]) override;
    void matTransVecMult(double x[], double r[]) override;
    void AddNewBlock(const int ci, const int pi, double * const data);
private:
    std::vector<MatrixBlock> m_Matrix;
};