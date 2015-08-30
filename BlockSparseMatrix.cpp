//
//  BlockSparseMatrix.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include <iostream>
#include "BlockSparseMatrix.h"

void BlockSparseMatrix::matVecMult(double x[], double r[])
{
    size_t i,j,k, n;
    for (k = 0, n = m_Matrix.size(); k < n; ++k)
    {
        //r[m_Matrix[k].pi * 2]       += *(m_Matrix[k].data) * x[m_Matrix[k].ci * 2];
        //r[(m_Matrix[k].pi * 2) + 1] += *(m_Matrix[k].data) * x[(m_Matrix[k].ci * 2) + 1];
        
        MatrixBlock block = (m_Matrix[k]);
        for(i=0;i<block.ilength;i++)
        {
            for(j=0;j<block.jlength;j++)
            {
                int i2 = block.ci+i;
                int j2 = block.pi+j;
                int idex = i*block.jlength+j;
                double val = block.data[idex];
                double prod = val * x[j2];
                r[i2]+=prod;
            }
        }
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])
{
    size_t i, n;
    for (i = 0, n = m_Matrix.size(); i < n; ++i)
    {
        r[m_Matrix[i].ci * 2] += *(m_Matrix[i].data) * x[m_Matrix[i].pi * 2];
        r[(m_Matrix[i].ci * 2) + 1] += *(m_Matrix[i].data) * x[(m_Matrix[i].pi * 2) + 1];
    }
}

void BlockSparseMatrix::AddNewBlock(const int ci, const int pi, const int ilength, const int jlength, const double * const data)
{
    m_Matrix.push_back(MatrixBlock(ci, pi, ilength, jlength, data));
}