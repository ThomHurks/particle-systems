#include "Solver.h"

Solver::Solver(const std::vector<Particle *> &pVector, const std::vector<Force *> &fVector,
               const std::vector<Force *> &cVector, double* CVector[],
               double* CDotVector[], BlockSparseMatrix &J, BlockSparseMatrix &JDot) :
m_ParticlesVector(pVector), m_ForcesVector(fVector), m_ConstraintsVector(cVector), m_CVector(CVector), m_CDotVector(CDotVector),
m_J(J), m_JDot(JDot)
{}

void Solver::simulation_step(const double dt, SolverType solverType)
{
    switch (solverType)
    {
        case SolverType::Euler:
            EulerSolver(dt);
            break;
        case SolverType::Midpoint:
            MidpointSolver(dt);
            break;
        case SolverType::RungeKutta4:
            RungeKutta4thOrderSolver(dt);
            break;
    }
}

void Solver::EulerSolver(const double dt)
{
    size_t ii, size = m_ParticlesVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    ParticleDerivative(dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position += dVector[ii].vec1; // XDot
        m_ParticlesVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
}

void Solver::MidpointSolver(const double dt)
{
    size_t ii, size = m_ParticlesVector.size();
    std::vector<Vec2fTuple> dVector(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    ParticleDerivative(dVector);
    ScaleVectorTuples(dVector, dt / 2); //notice the /2
    for (ii = 0; ii < size; ii++)
    {   //Make half an euler step, save original positions and velocities.
        originals[ii].vec1 = m_ParticlesVector[ii]->m_Position;
        originals[ii].vec2 = m_ParticlesVector[ii]->m_Velocity;
        m_ParticlesVector[ii]->m_Position += dVector[ii].vec1; // XDot
        m_ParticlesVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
    }
    ParticleDerivative(dVector);
    ScaleVectorTuples(dVector, dt);
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position = originals[ii].vec1 + dVector[ii].vec1; // XDot
        m_ParticlesVector[ii]->m_Velocity = originals[ii].vec2 + dVector[ii].vec2; // VDot
    }
}

void Solver::RungeKutta4thOrderSolver(const double dt)
{
    size_t ii, size = m_ParticlesVector.size();
    std::vector<Vec2fTuple> dVector1(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector2(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector3(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> dVector4(size); // Optimize this by reusing.
    std::vector<Vec2fTuple> originals(size);
    for (ii = 0; ii < size; ii++)
    {
        originals[ii].vec1 = m_ParticlesVector[ii]->m_Position;
        originals[ii].vec2 = m_ParticlesVector[ii]->m_Velocity;
    }
    //step 1
    ParticleDerivative(dVector1);
    ScaleVectorTuples(dVector1, dt);
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position = originals[ii].vec1 + dVector1[ii].vec1 /2; // XDot
        m_ParticlesVector[ii]->m_Velocity = originals[ii].vec2 + dVector1[ii].vec2 /2; // VDot
    }
    //step 2
    ParticleDerivative(dVector2);
    ScaleVectorTuples(dVector2, dt);
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position = originals[ii].vec1 + dVector2[ii].vec1 /2; // XDot
        m_ParticlesVector[ii]->m_Velocity = originals[ii].vec2 + dVector2[ii].vec2 /2; // VDot
    }
    //step 3
    ParticleDerivative(dVector3);
    ScaleVectorTuples(dVector3, dt);
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position = originals[ii].vec1 + dVector3[ii].vec1; // XDot
        m_ParticlesVector[ii]->m_Velocity = originals[ii].vec2 + dVector3[ii].vec2; // VDot
    }
    //step 4
    ParticleDerivative(dVector4);
    ScaleVectorTuples(dVector4, dt);
    //final positions & velocities
    for (ii = 0; ii < size; ii++)
    {
        m_ParticlesVector[ii]->m_Position = originals[ii].vec1 +
                (dVector1[ii].vec1 +2*dVector2[ii].vec1 +2*dVector3[ii].vec1 +dVector4[ii].vec1)/6; // XDot
        m_ParticlesVector[ii]->m_Velocity = originals[ii].vec2 +
                (dVector1[ii].vec2 +2*dVector2[ii].vec2 +2*dVector3[ii].vec2 +dVector4[ii].vec2)/6; // VDot
    }
    
}

void Solver::ParticleDerivative(std::vector<Vec2fTuple> & dVector)
{
    ClearForces(m_ParticlesVector);
    // First calculate forces:
    CalculateForces(m_ParticlesVector, m_ForcesVector);
    // Then calculate constraint forces:
    CalculateForces(m_ParticlesVector, m_ConstraintsVector);
    // Then derivatives:
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 = m_ParticlesVector[i]->m_Velocity; // XDot
        dVector[i].vec2 = m_ParticlesVector[i]->m_AccumulatedForce / m_ParticlesVector[i]->m_Mass; // VDot
    }
}

void Solver::Equation11(double ks, double kd)
{
    // Preparing variables for equation 11:
    size_t i;
    size_t n = m_ParticlesVector.size();

    // q contains all particle positions as 2-vectors. Size is n.
    Vec2f * q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { q[i] = m_ParticlesVector[i]->m_Position; }

    // qdot contains all velocities, but each dimension is stored separately as a double. The size is 2n.
    double * qdot = new double[2 * n];
    for (i = 0; i < n; ++i)
    {
        qdot[i * 2] = m_ParticlesVector[i]->m_Velocity[0];
        qdot[(i * 2) + 1] = m_ParticlesVector[i]->m_Velocity[1];
    }

    // M contains all particle masses as a double array of size n.
    double * M = new double[n];
    for (i = 0; i < n; ++i)
    { M[i] = m_ParticlesVector[i]->m_Mass; }

    // Q contains all accumulated forces as 2-vectors. Size is n.
    Vec2f * Q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { Q[i] = m_ParticlesVector[i]->m_AccumulatedForce; }

    // W contains the inverse of all particle masses as a double array of size n.
    double * W = new double[n];
    for (i = 0; i < n; ++i)
    { W[i] = 1 / m_ParticlesVector[i]->m_Mass; }

    // Start computing equation 11:

    // First calculate m_JDot times qdot. The result is a vector of size 2n since each dimension is stored separately.
    double * JDotqdot = new double[2 * n];
    std::fill(JDotqdot, JDotqdot + (2 * n), 0);
    m_JDot.matVecMult(qdot, JDotqdot);

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
    m_J.matVecMult(WQ, JWQ);

    // Then calculate ks times C. Size is n.
    double * ksC = new double[n];
    for (i = 0; i < n; ++i)
    { ksC[i] = ks * (*m_CVector)[i]; }

    // Then calculate kd times CDot. Size is n.
    double *kdCDot = new double[n];
    for (i = 0; i < n; ++i)
    { kdCDot[i] = kd * (*m_CDotVector)[i]; }

    // Now compute the entire right side of equation 11. The result is a double vector of size 2n.
    double * rightHandSide = new double[2 * n];
    for (i = 0; i < n; ++i)
    {
        rightHandSide[i * 2] = -JDotqdot[i] - JWQ[i] - ksC[i] - kdCDot[i];
        rightHandSide[(i * 2) + 1] = -JDotqdot[i + 1] - JWQ[i + 1] - ksC[i] - kdCDot[i];
    }

    // The left hand side of equation 11 is implemented inside the class JWJTranspose, an implicit matrix.
    double * lambda = new double[2 * n];
    JWJTranspose JWJTranspose(2 * n, W, m_J);

    int n_int = static_cast<int>((n * 2));
    int steps = 0; // 0 implies MAX_STEPS.
    ConjGrad(n_int, &JWJTranspose, rightHandSide, lambda, 0.1, &steps);
}

void Solver::ClearForces(const std::vector<Particle*> & pVector)
{
    size_t i, n;
    for (i = 0, n = pVector.size(); i < n; ++i)
    {
        pVector[i]->m_AccumulatedForce = Vec2f(0, 0);
    }
}

void Solver::CalculateForces(const std::vector<Particle*> & pVector, const std::vector<Force*> & fVector)
{
    size_t i, n;
    for (i = 0, n = fVector.size(); i < n; ++i)
    {
        fVector[i]->ApplyForce(pVector);
    }
}

void Solver::ScaleVectorTuples(std::vector<Vec2fTuple> &dVector, const double scaleFactor)
{
    size_t i, n;
    for (i = 0, n = dVector.size(); i < n; ++i)
    {
        dVector[i].vec1 *= scaleFactor;
        dVector[i].vec2 *= scaleFactor;
    }
}

