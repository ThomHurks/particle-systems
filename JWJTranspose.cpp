//
// Created by Thom Hurks on 20-08-15.
//

#include "JWJTranspose.h"

JWJTranspose::JWJTranspose(const size_t n, const double W[], BlockSparseMatrix &J) : n(n), W(W), J(J)
{
    m_IntermediateVector = new double[n];
}

void JWJTranspose::matVecMult(double x[], double r[])
{
    // First calculate JTranspose times X, store in m_IntermediateVector.
    J.matTransVecMult(x, m_IntermediateVector);
    // Then calculate W times JTransposeX, store in m_IntermediateVector.
    size_t i;
    for (i = 0; i < n; ++i)
    {
        m_IntermediateVector[i] *= W[i];
    }
    // Now do final multiplication with J to obtain JWJTransposeX.
    J.matVecMult(m_IntermediateVector, r);
}

void JWJTranspose::matTransVecMult(double x[], double r[])
{
    // Since the compound matrix JWJTranpose is symmetric.
    return matVecMult(x, r);
}

JWJTranspose::~JWJTranspose()
{
    delete[] m_IntermediateVector;
    m_IntermediateVector = nullptr;
}