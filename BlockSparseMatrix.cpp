//
//  BlockSparseMatrix.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include "BlockSparseMatrix.h"

void BlockSparseMatrix::matVecMult(double x[], double r[])
{
    int i;
    size_t n = m_Matrix.size();
    for (i = 0; i < n; ++i)
    {
        r[m_Matrix[i].pi] += *m_Matrix[i].data * x[m_Matrix[i].ci];
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])
{
    int i;
    size_t n = m_Matrix.size();
    for (i = 0; i < n; ++i)
    {
        r[m_Matrix[i].ci] += *m_Matrix[i].data * x[m_Matrix[i].pi];
    }
}

void BlockSparseMatrix::AddBlock(MatrixBlock block)
{
    m_Matrix.push_back(block);
}