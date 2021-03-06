#include <assert.h>
#include "Solver.h"
#include <stdlib.h>

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
        assert(!isnan(m_ParticlesVector[ii]->m_Position) && isfinite(m_ParticlesVector[ii]->m_Position));
        m_ParticlesVector[ii]->m_Velocity += dVector[ii].vec2; // VDot
        assert(!isnan(m_ParticlesVector[ii]->m_Velocity) && isfinite(m_ParticlesVector[ii]->m_Velocity));
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
    CalculateForces(m_ParticlesVector, m_ConstraintsVector);//What this does is fill J, Jdot, C and Cdot
    // Then solve for the constraint force on each particle: (this is equation 11)
    SolveConstraintForces(m_ks, m_kd, m_epsilon);
    // Then derivatives:
    size_t i, n;
    for (i = 0, n = dVector.size(); i < n; ++i)
    {
        Particle* p = m_ParticlesVector[i];
        dVector[i].vec1 = p->m_Velocity; // XDot
        dVector[i].vec2 = p->m_AccumulatedForce / p->m_Mass; // VDot
        assert(!isnan(dVector[i].vec1) && isfinite(dVector[i].vec1));
        assert(!isnan(dVector[i].vec2) && isfinite(dVector[i].vec2));
    }
}

// Equation 11
void Solver::SolveConstraintForces(const double ks, const double kd, const double epsilon)
{
    std::cout <<"J:"<<std::endl;
    m_J.print();
    std::cout << "Jdot"<<std::endl;
    m_JDot.print();
    
    // Preparing variables for equation 11:
    size_t i;
    const size_t n = m_ParticlesVector.size();
    const size_t two_n = 2 * n;
    const size_t m = m_ConstraintsVector.size();

    // qdot contains all velocities, but each dimension is stored separately as a double. The size is 2n.
    double *qdot = new double[two_n];
    for (i = 0; i < n; ++i)
    {
        qdot[i * 2] = m_ParticlesVector[i]->m_Velocity[0];
        qdot[(i * 2) + 1] = m_ParticlesVector[i]->m_Velocity[1];
    }
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(qdot[i]) && std::isfinite(qdot[i])); }

    // Create the C and CDot arrays that contain the constraint derivatives by gathering the C and CDot values.
    // The size of the arrays are n, since each constraint affects one or more particles.
    double *C = new double[m];
    double *CDot = new double[m];
    for (i = 0; i < m; ++i)
    {
        C[i] = (static_cast<Constraint*>(m_ConstraintsVector[i]))->GetC();
        CDot[i] = (static_cast<Constraint*>(m_ConstraintsVector[i]))->GetCDot();
    }
    for (i = 0; i < m; ++i) { assert(!std::isnan(C[i]) && std::isfinite(C[i])); }
    for (i = 0; i < m; ++i) { assert(!std::isnan(CDot[i]) && std::isfinite(CDot[i])); }

    // Q contains all accumulated forces as 2-vectors. Size is n.
    double *Q = new double[two_n];
    for (i = 0; i < n; ++i)
    { Q[2*i] = m_ParticlesVector[i]->m_AccumulatedForce[0];
    Q[2*i+1] = m_ParticlesVector[i]->m_AccumulatedForce[1];}
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(Q[i]) && std::isfinite(Q[i])); }

    // W contains the inverse of all particle masses as a double array of size n.
    double *W = new double[two_n];
    for (i = 0; i < n; ++i)
    { W[2*i] = 1.0 / m_ParticlesVector[i]->m_Mass;
    W[2*i+1] = 1.0 / m_ParticlesVector[i]->m_Mass;
    }
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(W[i]) && std::isfinite(W[i])); }

    // Start computing equation 11:

    // First calculate m_JDot times qdot. The result is a vector of size 2n since each dimension is stored separately.
    double *JDotqdot = new double[m];
    std::fill(JDotqdot, JDotqdot + m, 0.0);
    m_JDot.matVecMult(qdot, JDotqdot);
    for (i = 0; i < m; ++i) { assert(!std::isnan(JDotqdot[i]) && std::isfinite(JDotqdot[i])); }

    // Then calculate W times Q. The result is a vector of size 2n since each dimension is stored separately.
    double *WQ = new double[two_n];//W*Q = qDotDot
    for (i = 0; i < two_n; ++i)
    {
        WQ[i] = W[i] * Q[i];
    }
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(WQ[i]) && std::isfinite(WQ[i])); }

    // Then calculate J times WQ. The result is a vector of size m.
    double *JWQ = new double[m]; //= Cdotdot
    std::fill(JWQ, JWQ + m, 0.0);
    m_J.matVecMult(WQ, JWQ);
    for (i = 0; i < m; ++i) { assert(!std::isnan(JWQ[i]) && std::isfinite(JWQ[i])); }

    // Then calculate ks times C. C has size m
    double *ksC = new double[m];
    for (i = 0; i < m; ++i)
    {
        ksC[i] = ks * C[i];
    }
    for (i = 0; i < m; ++i) { assert(!std::isnan(ksC[i]) && std::isfinite(ksC[i])); }

    // Then calculate kd times CDot. CDot has size m.
    double *kdCDot = new double[m];
    for (i = 0; i < m; ++i)
    {
        kdCDot[i] = kd * CDot[i];
    }
    for (i = 0; i < m; ++i) { assert(!std::isnan(kdCDot[i]) && std::isfinite(kdCDot[i])); }

    // Now compute the entire right side of equation 11. The result is a double vector of size m.
    double *rightHandSide = new double[m];
    for (i = 0; i < m; ++i)
    {
        rightHandSide[i] = -JDotqdot[i] - JWQ[i] - ksC[i] - kdCDot[i];
    }
    
    printArray("JDotqdot",m,JDotqdot);
    printArray("JWQ",m,JWQ);
    printArray("ksC",m,ksC);
    printArray("kdCDot",m,kdCDot);
    printArray("RHS",m,rightHandSide);
    for (i = 0; i < m; ++i) { assert(!std::isnan(rightHandSide[i]) && std::isfinite(rightHandSide[i])); }

    // The left hand side of equation 11 is implemented inside the class JWJTranspose, an implicit matrix.
    double *lambda = new double[m];
    std::fill(lambda, lambda + m, 0.0);
    for (i = 0; i < m; ++i) { assert(!std::isnan(lambda[i]) && std::isfinite(lambda[i])); }
    JWJTranspose JWJTranspose(two_n, W, m_J);//two_n is the dimension of the intermediate vector, which is correct

    int m_int = static_cast<int>(m);
    int steps = 0; // 0 implies MAX_STEPS.
    std::cout << "Calling conjugate gradient algorithm...\n";
    double rSqrLen = ConjGrad(m_int, &JWJTranspose, lambda, rightHandSide, epsilon, &steps);// solve JWJT lambda = righthandside for lambda
    for (i = 0; i < m; ++i) { assert(!std::isnan(lambda[i]) && std::isfinite(lambda[i])); }
    std::cout << rSqrLen<<std::endl;
    std::cout << steps<<std::endl;
    double* JWJTlambda = new double[m];
    std::fill(JWJTlambda,JWJTlambda+m,0.0);
    JWJTranspose.matVecMult(lambda,JWJTlambda);
    for(i = 0; i < m; i++){JWJTlambda[i] = std::abs(JWJTlambda[i]-rightHandSide[i]);}
    printArray("JWJTl-RHS",m,JWJTlambda);
    double *QHat = new double[two_n]; //constraint forces -> length 2n
    std::fill(QHat, QHat + two_n, 0.0);
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(QHat[i]) && std::isfinite(QHat[i])); }
    m_J.matTransVecMult(lambda, QHat); //lambda has size m, so multiplying with JTrans yields a 2n vector
    for (i = 0; i < two_n; ++i) { assert(!std::isnan(QHat[i]) && std::isfinite(QHat[i])); }
    for(i = 0; i < n; ++i)
    {
        Vec2f constrainingForce = Vec2f(static_cast<float>(QHat[i * 2]), static_cast<float>(QHat[(i * 2) + 1]));
        assert(!isnan(constrainingForce) && isfinite(constrainingForce));
        m_ParticlesVector[i]->m_AccumulatedForce += constrainingForce;
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

void Solver::printArray(std::string name, int n, double* ar)
{
    std::cout<<name<<": ";
    int i;
    for(i = 0; i < n; i++)
    {
        std::cout<<ar[i]<<", ";
    }
    std::cout<<std::endl;
}

