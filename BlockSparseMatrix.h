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
    int ci; // constraint index
    int pi; // particle index
    double *data;
};

class BlockSparseMatrix : public implicitMatrixWithTrans
{
public:
    void matVecMult(double x[], double r[]);
    void matTransVecMult(double x[], double r[]);
    void AddBlock(MatrixBlock block);
private:
    std::vector<MatrixBlock> m_Matrix;
};