// TinkerToy.cpp : Defines the entry point for the console application.
//

#include "BlockSparseMatrix.h"
#include "Particle.h"
#include "Force.h"
#include "GravityForce.h"
#include "Solver.h"
#include "SpringForce.h"
#include "RodConstraint.h"
#include "CircularWireConstraint.h"
#include "imageio.h"
#include "MouseSpringForce.h"
#include "AngularSpring.h"
#include "FixedPointConstraint.h"

#include <vector>
#include <stdlib.h>
#include <stdio.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <iostream>

#else
#include <GL/glut.h>
#endif

/* global variables */

static int N;
static double dt, d;
static int dsim;
static int dump_frames;
static int frame_number;

static std::vector<Particle*> pVector;
static std::vector<Force*> fVector;
static std::vector<Force*> cVector;
static BlockSparseMatrix J, JDot;
static Solver * solver = nullptr;
static Solver::SolverType m_SolverType = Solver::SolverType::Euler;

static MouseSpringForce * msf = nullptr;
static bool userIsMouseInteracting = false;
static float * currentMousePosition = new float[2];

static int win_id;
static int win_x, win_y;
static int mouse_down[3];
static int mouse_release[3];
static int mouse_shiftclick[3];
static int omx, omy, mx, my;
static int hmx, hmy;

/*
----------------------------------------------------------------------
free/clear/allocate simulation data
----------------------------------------------------------------------
 */

static void free_data(void)
{
    // Clear() also calls destructors of the objects inside the vectors.
    size_t i, n;
    for (i = 0, n = pVector.size(); i < n; ++i)
    { delete pVector[i]; }
    pVector.clear();
    for (i = 0, n = fVector.size(); i < n; ++i)
    { delete fVector[i]; }
    fVector.clear();
    for (i = 0, n = cVector.size(); i < n; ++i)
    { delete cVector[i]; }
    cVector.clear();

    delete msf;
    msf = nullptr;

    // Make sure arrays are deleted with correct delete[] syntax:
    delete[] currentMousePosition;
    currentMousePosition = nullptr;

    delete solver;
    solver = nullptr;
}

static void init(void)
{
    /*
    BlockSparseMatrix M;
    double* *dataBlock1 = new double*[2];
    double* *dataBlock2 = new double*[2];
    double* *dataBlock3 = new double*[2];
    dataBlock1[0] = new double[1];
    dataBlock1[0][0]=1.0;
    dataBlock1[1] = new double[1];
    dataBlock1[1][0]=2.0;
    dataBlock2[0] = new double[1];
    dataBlock2[0][0]=3.0;
    dataBlock2[1] = new double[1];
    dataBlock2[1][0]=4.0;
    dataBlock3[0] = new double[1];
    dataBlock3[0][0]=5.0;
    dataBlock3[1] = new double[1];
    dataBlock3[1][0]=6.0;
    
    M.AddNewBlock(0,0,1,2,dataBlock1);
    M.AddNewBlock(1,0,1,2,dataBlock2);
    M.AddNewBlock(2,0,1,2,dataBlock3);
    std::cout<<"M:"<<std::endl;
    M.print();
    
    double* x = new double[3];
    x[0]=1.0;
    x[1]=2.0;
    x[2]=3.0;
    int i;
    
    std::cout<<"x: ";
    for(i=0; i < 3;i++)
    {
        std::cout<<x[i]<<" ";
    }
    std::cout<<std::endl;
    double* v = new double[2];
    std::fill(v,v+2,0);
    M.matVecMult(x,v);
    std::cout<<"v: ";
    for(i=0; i < 2;i++)
    {
        std::cout<<v[i]<<" ";
    }
    std::cout<<std::endl;
    
    
    double* x2 = new double[2];
    x2[0]=1.0;
    x2[1]=2.0;
    std::cout<<"x2: ";
    for(i=0; i < 2;i++)
    {
        std::cout<<x2[i]<<" ";
    }
    std::cout<<std::endl;
    double* v2 = new double[3];
    std::fill(v2,v2+3,0);
    M.matTransVecMult(x2,v2);
    std::cout<<"v2: ";
    for(i=0; i < 3;i++)
    {
        std::cout<<v2[i]<<" ";
    }
    std::cout<<std::endl;
    */
    
    
    
    
    // ks and kd are spring and damping constants for the constraint forces in equation 11.
    double ks = 1;
    double kd = 1;
    // epsilon is how low the linear conjugate gradient solver should go.
    double epsilon = 0.00000001;
    if (!solver)
    { solver = new Solver(pVector, fVector, cVector, J, JDot, ks, kd, epsilon); }

    if (!currentMousePosition)
    {
        currentMousePosition = new float[2];
    }
}

static void clear_data(void)
{
    size_t ii, size = pVector.size();

    for (ii = 0; ii < size; ii++) {
        pVector[ii]->reset();
    }
}

static void initTest(void)
{
    

    const float dist = 0.2;
    const Vec2f center(0.0, 0.0);
    const Vec2f offset(dist, 0.0);
    const Vec2f offset2(dist, dist);
    const Vec2f offset3(0, dist);

    // Create three particles, attach them to each other, then add a
    // circular wire constraint to the first.

    int particleID = 0;
    pVector.push_back(new Particle(center + offset, particleID++));
    pVector.push_back(new Particle(center + offset+ offset, particleID++));
    
    pVector.push_back(new Particle(center + offset+ offset+ offset, particleID++));

    // You should replace these with a vector generalized forces and one of
    // constraints...
    fVector.push_back(new GravityForce());
    fVector.push_back(new SpringForce(pVector[0],pVector[1],dist,1,1));
    int constraintID = 0;
    cVector.push_back(new CircularWireConstraint(pVector[0], center, dist, &J, &JDot, constraintID++));
    //cVector.push_back(new RodConstraint(pVector[0],pVector[1],dist,&J, &JDot, constraintID++));
    cVector.push_back(new FixedPointConstraint(pVector[2], center + offset+ offset+ offset, &J, &JDot, constraintID++));
    J.setDimensions(pVector.size(),cVector.size(),2);
    JDot.setDimensions(pVector.size(),cVector.size(),2);
    init();
}

static void initCloth(bool crossFibers)
{
    init();

    //Note; Without cross fibers appears to function better
    const float dist = 0.1f;
    const Vec2f topLeft(-0.75f, 0.75f);
    const Vec2f offset(dist, 0.0);
    const Vec2f offset2(0.0, -dist);
    const int dim = 15;

    int particleID = 0;
    int i;
    int j;
    //init particles
    for (i = 0; i <= dim; i++) {
        for (j = 0; j <= dim; j++) {
            pVector.push_back(new Particle(topLeft + offset * i + offset2*j, particleID++));
        }
    }
    double ks = 0.01;
    double kd = 0.01;
    double rest = 0.075;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            int cur = j * (dim + 1) + i;
            int right = cur + dim + 1;
            int below = cur + 1;
            fVector.push_back(new SpringForce(pVector[cur], pVector[right], rest, ks, kd));
            fVector.push_back(new SpringForce(pVector[cur], pVector[below], rest, ks, kd));
        }
    }
    for (i = 0; i < dim; i++) {
        int cur1 = (i + 1)*(dim + 1) - 1;
        int right = cur1 + dim + 1;
        fVector.push_back(new SpringForce(pVector[cur1], pVector[right], rest, ks, kd));

        int cur2 = i + dim * (dim + 1);
        int below = cur2 + 1;
        fVector.push_back(new SpringForce(pVector[cur2], pVector[below], rest, ks, kd));
    }
    if (crossFibers) {
        double drest = rest * sqrt(2);
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                int cur = j * (dim + 1) + i;
                int rightbelow = cur + dim + 2;
                fVector.push_back(new SpringForce(pVector[cur], pVector[rightbelow], drest, ks, kd));
            }
        }
        for (i = 1; i <= dim; i++) {
            for (j = 0; j < dim; j++) {
                int cur = j * (dim + 1) + i;
                int rightabove = cur + dim;
                fVector.push_back(new SpringForce(pVector[cur], pVector[rightabove], drest, ks, kd));
            }
        }
    }
}

static void initHair()
{
    init();

    const int internalParticles = 65; //amount of particles in hair is this + 2
    Vec2f start(0.25, -0.75f);
    Vec2f end(0.25, 0.0f);

    int particleID = 0;
    int i;
    Vec2f step = (end-start)/(internalParticles+1);
    for (i = 0; i <= internalParticles+1; i++) {
         pVector.push_back(new Particle(start + step * i, particleID++));
    }
    
    //m_ForcesVector.push_back(new GravityForce());
    float rest = magnitude(step);
    double ks = 0.05;
    double kd = 0.01;
    for (i = 0; i <= internalParticles; i++) {
         fVector.push_back(new SpringForce(pVector[i], pVector[i+1], rest, ks, kd));
    }
    ks = 0.01;
    kd = 0.1;
    double totalInternalAngle = 180.0f* (internalParticles);
    double angleDegrees = (totalInternalAngle/(internalParticles+2));
    double angleRadians = PI* angleDegrees/180.0f;
    for (i = 1; i <= internalParticles; i++) {
         fVector.push_back(new AngularSpring(pVector[i], pVector[i-1],pVector[i+1], angleRadians, ks, kd));
    }
    
}

/*
----------------------------------------------------------------------
OpenGL specific drawing routines
----------------------------------------------------------------------
 */

static void pre_display(void)
{
    glViewport(0, 0, win_x, win_y);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

static void post_display(void)
{
    // Write frames if necessary.
    if (dump_frames) {
        const int FRAME_INTERVAL = 4;
        if ((frame_number % FRAME_INTERVAL) == 0) {
            const int w = glutGet(GLUT_WINDOW_WIDTH);
            const int h = glutGet(GLUT_WINDOW_HEIGHT);
            if (w <= 0 || h <= 0) {
                exit(-1);
            }
            unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof (unsigned char));
            if (!buffer) {
                exit(-1);
            }
            // glRasterPos2i(0, 0);
            glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
            static char filename[80];
            sprintf(filename, "snapshots/img%.5i.png", frame_number / FRAME_INTERVAL);
            printf("Dumped %s.\n", filename);
            saveImageRGBA(filename, buffer, w, h);

            free(buffer);
        }
    }
    frame_number++;

    glutSwapBuffers();
}

static void draw_particles(void)
{
    size_t ii, size = pVector.size();

    for (ii = 0; ii < size; ii++) {
        pVector[ii]->draw();
    }
}

static void draw_forces(void)
{
    size_t i, n = fVector.size();
    for (i = 0; i < n; ++i) {
        fVector[i]->draw();
    }
}

static void draw_constraints(void)
{
    size_t i, n = cVector.size();
    for (i = 0; i < n; ++i) {
        cVector[i]->draw();
    }
}

/*
----------------------------------------------------------------------
relates mouse movements to tinker toy construction
----------------------------------------------------------------------
 */

static void get_from_UI()
{
    float i, j;
    // int size, flag;
    int hi, hj;
    // float x, y;
    if (!mouse_down[0] && !mouse_down[2] && !mouse_release[0]
            && !mouse_shiftclick[0] && !mouse_shiftclick[2]) return;

    i = ((mx / (float) win_x) * N);
    j = (((win_y - my) / (float) win_y) * N);

    if (i < 1 || i > N || j < 1 || j > N) return;


    hi = (int) ((hmx / (float) win_x) * N);
    hj = (int) (((win_y - hmy) / (float) win_y) * N);

    if (mouse_down[0]) {
        currentMousePosition[0] = 2 * (i - N / 2) / N;
        currentMousePosition[1] = 2 * (j - N / 2) / N;
        if (!userIsMouseInteracting) {
            userIsMouseInteracting = true;
            if (msf == nullptr) {
                msf = new MouseSpringForce(pVector[0], 0, 0.4, 0, currentMousePosition);
                fVector.push_back(msf);
            }
        }
    }

    if (mouse_down[2]) {
    }

    if (mouse_release[0]) {
        if (userIsMouseInteracting) {
            userIsMouseInteracting = false;
            if (fVector.size() > 0) {
                fVector.pop_back();
            }
            delete msf;
            msf = nullptr;
            currentMousePosition[0] = 0;
            currentMousePosition[1] = 0;
        }
        mouse_release[0] = 0; //only need to record this once
    }

    omx = mx;
    omy = my;
}

static void remap_GUI()
{
    size_t ii, size = pVector.size();
    for (ii = 0; ii < size; ii++) {
        pVector[ii]->reset();
    }
}

/*
----------------------------------------------------------------------
GLUT callback routines
----------------------------------------------------------------------
 */

static void key_func(unsigned char key, int x, int y)
{
    switch (key) {
        case 'c':
        case 'C':
            clear_data();
            break;

        case 'd':
        case 'D':
            dump_frames = !dump_frames;
            break;

        case 'q':
        case 'Q':
            free_data();
            exit(0); // implicit break

        case ' ':
            dsim = !dsim;
            break;
        case '1':
            m_SolverType = Solver::SolverType::Euler;
            break;
        case '2':
            m_SolverType = Solver::SolverType::Midpoint;
            break;
        case '3':
            m_SolverType = Solver::SolverType::RungeKutta4;
            break;
        case '!':
            dsim = false;
            free_data();
            initTest();
            break;
            case '@':
            dsim = false;
            free_data();
            initCloth(false);
            break;
        case '#':
            dsim = false;
            free_data();
            initCloth(true);
            break;
        case '$':
            dsim = false;
            free_data();
            initHair();
            break;
        default:
            std::cout << "Invalid input!\n";
            break;
    }
}

static void mouse_func(int button, int state, int x, int y)
{
    omx = mx = x;
    omx = my = y;

    if (!mouse_down[0]) {
        hmx = x;
        hmy = y;
    }
    if (mouse_down[button]) mouse_release[button] = state == GLUT_UP;
    if (mouse_down[button]) mouse_shiftclick[button] = glutGetModifiers() == GLUT_ACTIVE_SHIFT;
    mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func(int x, int y)
{
    mx = x;
    my = y;
}

static void reshape_func(int width, int height)
{
    glutSetWindow(win_id);
    glutReshapeWindow(width, height);

    win_x = width;
    win_y = height;
}

static void idle_func(void)
{
    if (dsim) {
        get_from_UI();
        solver->simulation_step(dt, m_SolverType);
    } else {
        get_from_UI();
        remap_GUI();
    }

    glutSetWindow(win_id);
    glutPostRedisplay();
}

static void display_func(void)
{
    pre_display();

    draw_forces();
    draw_constraints();
    draw_particles();

    post_display();
}

/*
----------------------------------------------------------------------
open_glut_window --- open a glut compatible window and set callbacks
----------------------------------------------------------------------
 */

static void open_glut_window(void)
{
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

    glutInitWindowPosition(0, 0);
    glutInitWindowSize(win_x, win_y);
    win_id = glutCreateWindow("Tinkertoys!");

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);

    pre_display();

    glutKeyboardFunc(key_func);
    glutMouseFunc(mouse_func);
    glutMotionFunc(motion_func);
    glutReshapeFunc(reshape_func);
    glutIdleFunc(idle_func);
    glutDisplayFunc(display_func);
}

/*
----------------------------------------------------------------------
main --- main routine
----------------------------------------------------------------------
 */

int main(int argc, char ** argv)
{
    glutInit(&argc, argv);

    if (argc == 1) {
        N = 64;
        dt = 0.1f;
        d = 5.f;
        fprintf(stderr, "Using defaults : N=%d dt=%g d=%g\n",
                N, dt, d);
    } else {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        d = atof(argv[3]);
    }

    printf("\n\nHow to use this application:\n\n");
    printf("\t Toggle construction/simulation display with the spacebar key\n");
    printf("\t Dump frames by pressing the 'd' key\n");
    printf("\t Change integratuion scheme by pressing 1,2  or 3 for Euler, Midpoint or RK4 respectively\n");
    printf("\t Change initial state by pressing !,@  or # \n");
    printf("\t Quit by pressing the 'q' key\n");

    dsim = 0;
    dump_frames = 0;
    frame_number = 0;

    initTest();

    win_x = 512;
    win_y = 512;
    open_glut_window();

    glutMainLoop();

    exit(0);
}

