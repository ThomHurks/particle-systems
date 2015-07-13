UNAME := $(shell uname)
CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = Solver.o Particle.o TinkerToy.o RodConstraint.o SpringForce.o CircularWireConstraint.o imageio.o

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
