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
    for (k = 0, n = m_Matrix.size(); k < n; ++k) {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i) {
            size_t gi = block.ci + i;
            for (j = 0; j < block.jlength; ++j)//2
            {
                // Todo: ensure correct indexing into block. done
                size_t gj = block.pj * 2 + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gj];
                r[gi] += prod;
            }
        }
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])//x.length is equal to the width (greatest ci+ilength) of the maatrix
{
    size_t i, j, k, n;
    for (k = 0, n = m_Matrix.size(); k < n; ++k) {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i) {
            size_t gi = block.ci + i;
            for (j = 0; j < block.jlength; ++j) {
                // Todo: ensure correct indexing into block.
                size_t gj = 2 * block.pj + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gi];
                r[gj] += prod;
            }
        }
    }
}

void BlockSparseMatrix::AddNewBlock(const int ci, const int pi, const int ilength, const int jlength, double* const data[])
{
    m_Matrix.push_back(MatrixBlock(ci, pi, ilength, jlength, data));
}

void BlockSparseMatrix::print()
{
    int maxci = 0;
    int maxpj = 0;
    int k;
    for (k = 0; k < m_Matrix.size(); ++k) {
        maxci = std::max(maxci, m_Matrix[k].ci);
        maxpj = std::max(maxpj, m_Matrix[k].pj);
    }
    int dimi =maxci+1;
    int dimj = 2*(maxpj+1);
    double board[dimi][dimj]; //creates a 9*9 matrix or a 2d array.

    for (int i = 0; i < dimi; i++) //This loops on the rows.
    {
        for (int j = 0; j < dimj; j++) //This loops on the columns
        {
            board[i][j] = 0; //replace by some fill command
        }
    }
    int i, j;
    for (k = 0; k < m_Matrix.size(); ++k) {
        for (i = 0; i < m_Matrix[k].ilength; i++) {
            int gi = m_Matrix[k].ci+i;
            for (j = 0; j < m_Matrix[k].jlength; j++) {
                int gj = m_Matrix[k].pj*2+j;
                board[gi][gj]=*(m_Matrix[k].data[m_Matrix[k].Index(i,j)]);
            }
        }
    }
    for(int j=0; j<dimj; j++) 
	{
		for(int i=0; i<dimi; i++)
		{
			std::cout << board[i][j]  << "  ";
		}
		std::cout << std::endl;
	}
}
