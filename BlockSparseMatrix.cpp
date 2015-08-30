//
//  BlockSparseMatrix.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include <iostream>
#include "BlockSparseMatrix.h"

void BlockSparseMatrix::matVecMult(double x[], double r[])//x.length is equal (or greater) than the heigh (greatest pi+jlength) of the maatrix
{
    size_t i, j, k, n;
    for (k = 0, n = m_Matrix.size(); k < n; ++k)
    {
        // Todo: fix this as data is no longer a pointer to a double, but a matrix of pointers to doubles.
        //r[m_Matrix[k].pi * 2]       += *(m_Matrix[k].data) * x[m_Matrix[k].ci * 2];
        //r[(m_Matrix[k].pi * 2) + 1] += *(m_Matrix[k].data) * x[(m_Matrix[k].ci * 2) + 1];
        
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i)
        {
            size_t gi = block.ci + i;
            for (j = 0; j < block.jlength; ++j)
            {
                // Todo: ensure correct indexing into block. done
                size_t gj = block.pi + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gi];
                r[gj]+=prod;
            }
        }
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])//x.length is equal (or greater) than the width (greatest ci+ilength) of the maatrix
{
    size_t i,j,k, n;
    for (k = 0, n = m_Matrix.size(); k < n; ++i)
    {
        // Todo: fix this as data is no longer a pointer to a double, but a matrix of pointers to doubles.
        //r[m_Matrix[i].ci * 2] += *(m_Matrix[i].data) * x[m_Matrix[i].pi * 2];
        //r[(m_Matrix[i].ci * 2) + 1] += *(m_Matrix[i].data) * x[(m_Matrix[i].pi * 2) + 1];
        
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i)
        {
            size_t gi = block.ci + i;
            for (j = 0; j < block.jlength; ++j)
            {
                // Todo: ensure correct indexing into block.
                size_t gj = block.pi + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gj];
                r[gi]+=prod;
            }
        }
    }
}

void BlockSparseMatrix::AddNewBlock(const int ci, const int pi, const int ilength, const int jlength, double* const data[])
{
    m_Matrix.push_back(MatrixBlock(ci, pi, ilength, jlength, data));
}