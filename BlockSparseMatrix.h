//
//  BlockSparseMatrix.h
//  Particle Systems
//
//  Created by Thom Hurks on 28-07-15.
//  Copyright (c) 2015 Thom Hurks. All rights reserved.
//

#include <vector>

struct MatrixBlock {
    int i;
    int j;
    int ilength;
    int jlength;
    float *data;
};

class BlockSparseMatrix : public implicitMatrixWithTrans
{
public:
    void matVecMult(double x[], double r[]);
    void matTransVecMult(double x[], double r[]);
private:
    std::vector<MatrixBlock> m_Matrix;
};