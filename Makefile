

CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H 
OBJS = Solver.o Particle.o TinkerToy.o RodConstraint.o SpringForce.o CircularWireConstraint.o imageio.o

project1: $(OBJS)
	$(CXX) -o $@ $^ -lpng -framework OpenGL -framework GLUT
clean:
	rm $(OBJS) project1
