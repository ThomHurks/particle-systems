#include "Solver.h"

void simulation_step(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                     const std::vector<Force*> & cVector, const SolverType m_SolverType, double* CVector[],
                     double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt)
{
    switch (m_SolverType)
    {
        case SolverType::Euler:
            EulerSolver(pVector, fVector, cVector, CVector, CDotVector, J, JDot,  dt);
            break;
        case SolverType::Midpoint:
            MidpointSolver(pVector, fVector, cVector, CVector, CDotVector, J, JDot, dt);
            break;
        case SolverType::RungeKutta4:
            RungeKutta4thOrderSolver(pVector, fVector, cVector, CVector, CDotVector, J, JDot, dt);
            break;
    }
}

void EulerSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                 const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                 BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position += dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
}

void MidpointSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                    const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                    BlockSparseMatrix * J, BlockSparseMatrix * JDot, const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt / 2); //notice the /2
    for (ii = 0; ii < size; ii++)
    {   //Make half an euler step, save original positions and velocities.
        originals[ii].vec1 = pVector[ii]->m_Position;
        originals[ii].vec2 = pVector[ii]->m_Velocity;
        pVector[ii]->m_Position += dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
    ParticleDerivative(pVector, fVector, cVector, dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector[ii].vec1; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector[ii].vec2; // VDot
    }
}

void RungeKutta4thOrderSolver(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                              const std::vector<Force*> & cVector, double* CVector[], double* CDotVector[],
                              BlockSparseMatrix * J, BlockSparseMatrix * JDot,  const double dt)
{
    size_t ii, size = pVector.size();
    std::vector<Vec2fTuple> dVector1(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector2(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector3(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector4(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    for (ii = 0; ii < size; ii++)
    {
        originals[ii].vec1 = pVector[ii]->m_Position;
        originals[ii].vec2 = pVector[ii]->m_Velocity;
    }
    //step 1
    ParticleDerivative(pVector, fVector, cVector, dVector1);
    ScaleVectorTuples(dVector1, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector1[ii].vec1 /2; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector1[ii].vec2 /2; // VDot
    }
    //step 2
    ParticleDerivative(pVector, fVector, cVector, dVector2);
    ScaleVectorTuples(dVector2, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector2[ii].vec1 /2; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector2[ii].vec2 /2; // VDot
    }
    //step 3
    ParticleDerivative(pVector, fVector, cVector, dVector3);
    ScaleVectorTuples(dVector3, dt);
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 + dVector3[ii].vec1; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 + dVector3[ii].vec2; // VDot
    }
    //step 4
    ParticleDerivative(pVector, fVector, cVector, dVector4);
    ScaleVectorTuples(dVector4, dt);
    //final positions & velocities
    for (ii = 0; ii < size; ii++)
    {
        pVector[ii]->m_Position = originals[ii].vec1 +
                (dVector1[ii].vec1 +2*dVector2[ii].vec1 +2*dVector3[ii].vec1 +dVector4[ii].vec1)/6; // XDot
        pVector[ii]->m_Velocity = originals[ii].vec2 +
                (dVector1[ii].vec2 +2*dVector2[ii].vec2 +2*dVector3[ii].vec2 +dVector4[ii].vec2)/6; // VDot
    }
    
}

void ParticleDerivative(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector,
                        const std::vector<Force*> & cVector, std::vector<Vec2fTuple> & dVector)
{
    ClearForces(pVector);
    // First calculate forces:
    CalculateForces(pVector, fVector);
    // Then calculate constraint forces:
    CalculateForces(pVector, cVector);
    // Then derivatives:
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 = pVector[i]->m_Velocity; // XDot
        dVector[i].vec2 = pVector[i]->m_AccumulatedForce / pVector[i]->m_Mass; // VDot
    }
}

void Equation11(const std::vector<Particle*> & pVector, const std::vector<Force*> & cVector, double* CVector[],
                double* CDotVector[], BlockSparseMatrix * J, BlockSparseMatrix * JDot)
{
    // Todo: pass in ks and kd as parameters.
    double ks, kd;
    ks = kd = 1;

    // Preparing variables for equation 11:
    size_t i, n = pVector.size();

    // q contains all particle positions as 2-vectors. Size is n.
    Vec2f * q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { q[i] = pVector[i]->m_Position; }

    // qdot contains all velocities, but each dimension is stored separately as a double. The size is 2n.
    double * qdot = new double[2 * n];
    for (i = 0; i < n; ++i)
    {
        qdot[i * 2] = pVector[i]->m_Velocity[0];
        qdot[(i * 2) + 1] = pVector[i]->m_Velocity[1];
    }

    // M contains all particle masses as a double array of size n.
    double * M = new double[n];
    for (i = 0; i < n; ++i)
    { M[i] = pVector[i]->m_Mass; }

    // Q contains all accumulated forces as 2-vectors. Size is n.
    Vec2f * Q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { Q[i] = pVector[i]->m_AccumulatedForce; }

    // W contains the inverse of all particle masses as a double array of size n.
    double * W = new double[n];
    for (i = 0; i < n; ++i)
    { W[i] = 1 / pVector[i]->m_Mass; }

    // Start computing equation 11:

    // First calculate JDot times qdot. The result is a vector of size 2n since each dimension is stored separately.
    double * JDotqdot = new double[2 * n];
    std::fill(JDotqdot, JDotqdot + (2 * n), 0);
    JDot->matVecMult(qdot, JDotqdot);

    // Then calculate W times Q. The result is a vector of size 2n since each dimension is stored separately.
    double * WQ = new double[2 * n];
    for (i = 0; i < n; ++i)
    {
        WQ[i * 2] = W[i] * Q[i][0];
        WQ[(i * 2) + 1] = W[i] * Q[i][1];
    }

    // Then calculate J times WQ. The result is a vector of size 2n since each dimension is stored separately.
    double * JWQ = new double[2 * n];
    std::fill(JWQ, JWQ + (2 * n), 0);
    J->matVecMult(WQ, JWQ);

    // Then calculate ks times C. Size is n.
    double * ksC = new double[n];
    for (int i = 0; i < n; ++i)
    { ksC[i] = ks * (*CVector)[i]; }

    // Then calculate kd times CDot. Size is n.
    double *kdCDot = new double[n];
    for (int i = 0; i < n; ++i)
    { kdCDot[i] = kd * (*CDotVector)[i]; }

    // Now compute the entire right side of equation 11. The result is a double vector of size 2n.
    double * rightHandSide = new double[2 * n];
    for (int i = 0; i < n; ++i)
    {
        rightHandSide[i * 2] = -JDotqdot[i] - JWQ[i] - ksC[i] - kdCDot[i];
        rightHandSide[(i * 2) + 1] = -JDotqdot[i + 1] - JWQ[i + 1] - ksC[i] - kdCDot[i];
    }

    // I have no idea how to get the left side of equation 11. The course notes do not explain it clearly.
}

void ClearForces(const std::vector<Particle*> & pVector)
{
    size_t i, n = pVector.size();
    for (i = 0; i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce = Vec2f(0, 0);
    }
}

void CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector)
{
    size_t i, n = fVector.size();
    for (i = 0; i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

void ScaleVectorTuples(std::vector<Vec2fTuple> &dVector, const double scaleFactor)
{
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 *= scaleFactor;
        dVector[i].vec2 *= scaleFactor;
    }
}

