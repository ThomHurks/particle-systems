#include "Solver.h"

Solver::Solver(const std::vector<Particle *> &pVector, const std::vector<Force *> &fVector,
               const std::vector<Force *> &cVector, BlockSparseMatrix &J, BlockSparseMatrix &JDot, const double ks,
               const double kd, const double epsilon) :
m_ParticlesVector(pVector), m_ForcesVector(fVector), m_ConstraintsVector(cVector), m_J(J), m_JDot(JDot), m_ks(ks),
m_kd(kd), m_epsilon(epsilon)
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
    // Then solve for the constraint force on each particle: (this is equation 11)
    SolveConstraintForces(m_ks, m_kd, m_epsilon);
    // Then derivatives:
    size_t i, n = dVector.size();
    for (i = 0; i < n; ++i)
    {
        dVector[i].vec1 = m_ParticlesVector[i]->m_Velocity; // XDot
        dVector[i].vec2 = m_ParticlesVector[i]->m_AccumulatedForce / m_ParticlesVector[i]->m_Mass; // VDot
    }
}

// Equation 11
void Solver::SolveConstraintForces(const double ks, const double kd, const double epsilon)
{
    // Preparing variables for equation 11:
    size_t i;
    size_t n = m_ParticlesVector.size();
    size_t two_n = 2 * n;

    // qdot contains all velocities, but each dimension is stored separately as a double. The size is 2n.
    double *qdot = new double[2 * n];
    for (i = 0; i < n; ++i)
    {
        qdot[i * 2] = m_ParticlesVector[i]->m_Velocity[0];
        qdot[(i * 2) + 1] = m_ParticlesVector[i]->m_Velocity[1];
    }

    // Create the C and CDot arrays that contain the constraint derivatives by gathering the C and CDot values.
    // The size of the arrays are n, since each constraint affects one or more particles.
    size_t n_Constraints = m_ConstraintsVector.size();
    double *C = new double[n_Constraints];
    double *CDot = new double[n_Constraints];
    for (i = 0; i < n_Constraints; ++i)
    {
        C[i] = (static_cast<Constraint*>(m_ConstraintsVector[i]))->GetC();
        CDot[i] = (static_cast<Constraint*>(m_ConstraintsVector[i]))->GetCDot();
    }

    // Q contains all accumulated forces as 2-vectors. Size is n.
    Vec2f *Q = new Vec2f[n];
    for (i = 0; i < n; ++i)
    { Q[i] = m_ParticlesVector[i]->m_AccumulatedForce; }

    // W contains the inverse of all particle masses as a double array of size n.
    double *W = new double[n];
    for (i = 0; i < n; ++i)
    { W[i] = 1 / m_ParticlesVector[i]->m_Mass; }

    // Start computing equation 11:

    // First calculate m_JDot times qdot. The result is a vector of size 2n since each dimension is stored separately.
    double *JDotqdot = new double[two_n];
    std::fill(JDotqdot, JDotqdot + two_n, 0);
    m_JDot.matVecMult(qdot, JDotqdot);

    // Then calculate W times Q. The result is a vector of size 2n since each dimension is stored separately.
    double *WQ = new double[two_n];
    for (i = 0; i < n; ++i)
    {
        WQ[i * 2] = W[i] * Q[i][0];
        WQ[(i * 2) + 1] = W[i] * Q[i][1];
    }

    // Then calculate J times WQ. The result is a vector of size 2n since each dimension is stored separately.
    double *JWQ = new double[two_n];
    std::fill(JWQ, JWQ + two_n, 0);
    m_J.matVecMult(WQ, JWQ);

    // Then calculate ks times C. C has size n, but we need 2n, so duplicate all values.
    double *ksC = new double[two_n];
    for (i = 0; i < n; ++i)
    {
        ksC[i * 2] = ks * C[i];
        ksC[(i * 2) + 1] = ks * C[i];
    }

    // Then calculate kd times CDot. CDot has size n, but we need 2n, so duplicate all values.
    double *kdCDot = new double[two_n];
    for (i = 0; i < n; ++i)
    {
        kdCDot[i * 2] = kd * CDot[i];
        kdCDot[(i * 2) + 1] = kd * CDot[i];
    }

    // Now compute the entire right side of equation 11. The result is a double vector of size 2n.
    double *rightHandSide = new double[two_n];
    for (i = 0; i < two_n; ++i)
    {
        rightHandSide[i] = -JDotqdot[i] - JWQ[i] - ksC[i] - kdCDot[i];
    }

    // The left hand side of equation 11 is implemented inside the class JWJTranspose, an implicit matrix.
    double *lambda = new double[two_n];
    std::fill(lambda, lambda + two_n, 0);
    JWJTranspose JWJTranspose(two_n, W, m_J);

    int n_int = static_cast<int>(two_n);
    int steps = 0; // 0 implies MAX_STEPS.
    std::cout << "Calling conjugate gradient algorithm...\n";
    double rSqrLen = ConjGrad(n_int, &JWJTranspose, rightHandSide, lambda, epsilon, &steps);
    std::cout << rSqrLen;
    std::cout << '\n';
    std::cout << steps;
    std::cout << '\n';

    double *QHat = new double[two_n];
    std::fill(QHat, QHat + two_n, 0);
    m_J.matTransVecMult(lambda, QHat);
    for(i = 0; i < n; ++i)
    {
        Vec2f constrainingForce = Vec2f(static_cast<float>(QHat[i * 2]), static_cast<float>(QHat[(i * 2) + 1]));
        if (!isnan(constrainingForce))
        { m_ParticlesVector[i]->m_AccumulatedForce += constrainingForce; }
    }

    // Free all memory.
    delete[] qdot;
    delete[] C;
    delete[] CDot;
    delete[] Q;
    delete[] W;
    delete[] JDotqdot;
    delete[] WQ;
    delete[] JWQ;
    delete[] ksC;
    delete[] kdCDot;
    delete[] rightHandSide;
    delete[] lambda;
    delete[] QHat;
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

