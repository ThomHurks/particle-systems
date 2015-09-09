//
//  BlockSparseMatrix.cpp
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include <iostream>
#include <assert.h>
#include "BlockSparseMatrix.h"

void BlockSparseMatrix::setDimensions(const int n, const int m, const int d)
{
    
    m_n=n;
    m_m=m;
    m_d=d;
    m_cWidths = new int[m];
    int i;
    for(i=0;i<m_Matrix.size();i++)
    {
        m_cWidths[m_Matrix[i].ci]=m_Matrix[i].ilength;
    }
    m_pHeights = new int[n];
    std::fill(m_pHeights,m_pHeights+n,d);
    m_gJs = new int[m];
    m_gIs = new int[n];
    if(n>0){m_gJs[0]=0;}
    if(m>0){m_gIs[0]=0;}
    for(i = 1; i<m;i++)
    {
        m_gIs[i]=m_gIs[i-1]+m_cWidths[i-1];
    }
    for(i = 1; i < n; i++)
    {
        m_gJs[i]=m_gJs[i-1]+m_pHeights[i-1];
    }
    m_cWidth = m_gIs[m-1]+m_cWidths[m-1];
    m_pHeight = m_gJs[n-1]+m_pHeights[n-1];
}

int BlockSparseMatrix::getGlobalI(const int i)
{
    return m_gIs[i];
}


int BlockSparseMatrix::getGlobalJ(const int j)
{
    return m_gJs[j];
}


void BlockSparseMatrix::matVecMult(double x[], double r[]) // x.length is equal (or greater) than the height (greatest pi+jlength) of the matrix
{
    size_t i, j, k, n;
    //r should be m 0's, x should be n doubles
    for(i = 0; i < m_cWidth;i++){r[i]=0.0;}
    for(i = 0; i < m_pHeight;i++){assert(!isnan(x[i]) && isfinite(x[i]));}
    for (k = 0, n = m_Matrix.size(); k < n; ++k) {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i) {
            size_t gi = getGlobalI(block.ci) + i;
            for (j = 0; j < block.jlength; ++j)//2
            {
                // Todo: ensure correct indexing into block. done
                size_t gj = getGlobalJ(block.pj) + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gj];
                assert(!isnan(prod) && isfinite(prod));
                r[gi] += prod;
                assert(!isnan(r[gi]) && isfinite(r[gi]));
            }
        }
    }
}

void BlockSparseMatrix::matTransVecMult(double x[], double r[])//x.length is equal to the width (greatest ci+ilength) of the maatrix
{
    size_t i, j, k, n;
    
    //r should be m 0's, x should be n doubles
    for(i = 0; i < m_pHeight;i++){r[i]=0.0;}
    for(i = 0; i < m_cWidth;i++){assert(!isnan(x[i]) && isfinite(x[i]));}
    for (k = 0, n = m_Matrix.size(); k < n; ++k) {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i) {
            size_t gi = getGlobalI(block.ci) + i;
            for (j = 0; j < block.jlength; ++j) {
                // Todo: ensure correct indexing into block.
                size_t gj = getGlobalJ(block.pj) + j;
                double val = *(block.data[block.Index(i, j)]);
                double prod = val * x[gi];
                assert(!isnan(prod) && isfinite(prod));
                r[gj] += prod;
                assert(!isnan(r[gj]) && isfinite(r[gj]));
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
    std::cout<<"WxH = " << m_cWidth <<"x"<<m_pHeight<<std::endl;
    for(int j=0; j<dimj; j++) 
	{
		for(int i=0; i<dimi; i++)
		{
			std::cout << board[i][j]  << "  ";
		}
		std::cout << std::endl;
	}
}
