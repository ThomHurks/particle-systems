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
    const int ci;       // constraint index (i)
    const int pj;       // particle index (j)
    const int ilength;  // width of block (always 1 in our case), may vary from constraint to constraint
    const int jlength;  // height of block (always 2 in our case) = dimension of particle space
    double* const *data;// array of pointers to doubles

    MatrixBlock(const int ci, const int pj, const int ilength, const int jlength, double* const data[]) :
    ci(ci), pj(pj), ilength(ilength), jlength(jlength), data(data) {}

    size_t Index(const size_t i, const size_t j) const { return (j * ilength) + i;}
    
};

class BlockSparseMatrix : public implicitMatrixWithTrans
{
public:
    void setDimensions(const int n, const int m, const int d);
    void matVecMult(double x[], double r[]) override;
    void matTransVecMult(double x[], double r[]) override;
    void AddNewBlock(const int ci, const int pj, const int ilength, const int jlength, double* const data[]);
    void print();
    int getGlobalI(const int i);
    int getGlobalJ(const int j);
    void empty();
private:
    std::vector<MatrixBlock> m_Matrix;
    int m_n;//the ammount of particles in the system
    int m_m;//the ammount of constraints in the system
    int m_d;//the dimension of the system
    int *m_cWidths;//at index i, the width of constraint i is stored
    int *m_pHeights;//at index j, the height of particle j (d) is stored
    int m_cWidth;//the width of the total matrix
    int m_pHeight;//the height of the total matrix
    int *m_gIs;//the combined width of the constraints up to i-1
    int *m_gJs;//the combined height of the particles up to j-1
    bool m_initialised = false;
};