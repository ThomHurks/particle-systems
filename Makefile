UNAME := $(shell uname)
CXX = g++
CXXFLAGS = -g -Wall -Wextra -Wmissing-declarations -Wnull-dereference -Wswitch-default -Wswitch-enum -Wuninitialized -Wbad-function-cast -Wabstract-final-class -Wfloat-equal -Wcast-qual -Woverloaded-virtual -Wno-sign-compare -Wno-deprecated-declarations -Iinclude -DHAVE_CONFIG_H -std=c++0x

OBJS = Solver.o Particle.o TinkerToy.o RodConstraint.o SpringForce.o GravityForce.o CircularWireConstraint.o imageio.o BlockSparseMatrix.o AngularSpring.o MouseSpringForce.o JWJTranspose.o linearSolver.o

ifeq ($(UNAME), Darwin)
	TARGET = project1
else
	TARGET = project1.exe
endif

project1: $(OBJS)
	$(CXX) -o $@ $^ -lpng -framework OpenGL -framework GLUT
project1.exe: $(OBJS)
	$(CXX) -o $@ $^ -lpng -lGL -lGLU -lglut
clean:
	rm $(OBJS) project1
