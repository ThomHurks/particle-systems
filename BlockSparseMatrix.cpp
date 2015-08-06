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
        r[m_Matrix[i].pi * 2] += *m_Matrix[i].data * x[m_Matrix[i].ci * 2];
        r[(m_Matrix[i].pi * 2) + 1] += *m_Matrix[i].data * x[(m_Matrix[i].ci * 2) + 1];
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])
{
    int i;
    size_t n = m_Matrix.size();
    for (i = 0; i < n; ++i)
    {
        r[m_Matrix[i].ci * 2] += *m_Matrix[i].data * x[m_Matrix[i].pi * 2];
        r[(m_Matrix[i].ci * 2) + 1] += *m_Matrix[i].data * x[(m_Matrix[i].pi * 2) + 1];
    }
}

void BlockSparseMatrix::AddNewBlock(const int ci, const int pi, double * const data)
{
    m_Matrix.push_back(MatrixBlock(ci, pi, data));
}