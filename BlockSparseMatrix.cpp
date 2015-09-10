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

BlockSparseMatrix::BlockSparseMatrix(const size_t numParticles, const size_t numConstraints, const int particleDimension)
: m_ParticleDimension(particleDimension)
{
    SetDimensions(numParticles, numConstraints);
    if (numConstraints <= 0 || numParticles <= 0)
    { m_initialised = false; } // Need to explicitly call SetDimensions first.
}

BlockSparseMatrix::~BlockSparseMatrix()
{
    Clear();
}

void BlockSparseMatrix::Clear()
{
    m_Matrix.clear();
    delete[] m_gJs;
    m_gJs = nullptr;
    delete[] m_gIs;
    m_gIs = nullptr;
    m_cWidth = 0;
    m_pHeight = 0;
    m_initialised = false;
}

void BlockSparseMatrix::SetDimensions(const size_t numParticles, const size_t numConstraints)
{
    int *m_cWidths = new int[numConstraints];   // at index i, the width of constraint i is stored
    int *m_pHeights = new int[numParticles];    // at index j, the height of particle j (d) is stored

    std::fill(m_cWidths, m_cWidths + numConstraints, 0);
    std::fill(m_pHeights, m_pHeights + numParticles, m_ParticleDimension);

    size_t i, k;
    // Fixme: what if a constraint has multiple matrix blocks? This can happen according to course notes.
    // Fixme: The for loop will then overwrite some values as the blocks will have the same constraint index.
    assert(m_Matrix.size() >= numConstraints); // Each constraint creates at least one matrix block.
    for (i = 0, k = m_Matrix.size(); i < k; ++i)
    {
        assert(m_Matrix[i].ci < numConstraints);
        m_cWidths[m_Matrix[i].ci] = m_Matrix[i].ilength;
    }

    // Initialize the "base cases":
    delete[] m_gJs;
    m_gJs = new int[numParticles];
    if (numParticles > 0)   { m_gJs[0] = 0; }
    delete[] m_gIs;
    m_gIs = new int[numConstraints];
    if (numConstraints > 0) { m_gIs[0] = 0; }

    for (i = 1; i < numConstraints; ++i)    { m_gIs[i] = m_gIs[i - 1] + m_cWidths[i - 1];  }
    for (i = 1; i < numParticles; ++i)      { m_gJs[i] = m_gJs[i - 1] + m_pHeights[i - 1]; }

    if (numConstraints > 0)
    { m_cWidth =  m_gIs[numConstraints - 1] + m_cWidths[numConstraints - 1]; }
    else
    { m_cWidth = 0; }
    if (numParticles > 0)
    { m_pHeight = m_gJs[numParticles - 1]   + m_pHeights[numParticles - 1]; } // Todo: Isn't pHeights always the dimension?
    else
    { m_pHeight = 0; }
    m_initialised = true;

    delete[] m_cWidths;
    delete[] m_pHeights;
}

void BlockSparseMatrix::matVecMult(double x[], double r[]) // x.length is equal (or greater) than the height (greatest pi+jlength) of the matrix
{
    assert(m_initialised);
    size_t i, j, k, n;
    // r should be m 0's, x should be n doubles
    for (i = 0; i < m_cWidth; i++)
    {
        r[i] = 0.0;
    }
    for (i = 0; i < m_pHeight; i++)
    {
        assert(!isnan(x[i]) && isfinite(x[i]));
    }
    for (k = 0, n = m_Matrix.size(); k < n; ++k)
    {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i)
        {
            size_t gi = getGlobalI(block.ci) + i;
            for (j = 0; j < block.jlength; ++j) // 2
            {
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
    assert(m_initialised);
    size_t i, j, k, n;
    // r should be m 0's, x should be n doubles
    for (i = 0; i < m_pHeight; i++)
    {
        r[i] = 0.0;
    }
    for (i = 0; i < m_cWidth; i++)
    {
        assert(!isnan(x[i]) && isfinite(x[i]));
    }
    for (k = 0, n = m_Matrix.size(); k < n; ++k)
    {
        MatrixBlock block = m_Matrix[k];
        for (i = 0; i < block.ilength; ++i)
        {
            size_t gi = getGlobalI(block.ci) + i;
            for (j = 0; j < block.jlength; ++j)
            {
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

void BlockSparseMatrix::print() const
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
