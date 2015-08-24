cmake_minimum_required(VERSION 3.2)
project(tue_scg_particle_systems)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -Wall -Wno-sign-compare -Woverloaded-virtual -Iinclude -DHAVE_CONFIG_H")

set(SOURCE_FILES
        include/gfx/config.h
        include/gfx/geom3d.h
        include/gfx/geom4d.h
        include/gfx/gfx.h
        include/gfx/mat2.h
        include/gfx/mat3.h
        include/gfx/mat4.h
        include/gfx/quat.h
        include/gfx/symmat2.h
        include/gfx/symmat3.h
        include/gfx/symmat4.h
        include/gfx/vec2.h
        include/gfx/vec3.h
        include/gfx/vec4.h
        BlockSparseMatrix.cpp
        BlockSparseMatrix.h
        CircularWireConstraint.cpp
        CircularWireConstraint.h
        Force.h
        GravityForce.cpp
        GravityForce.h
        imageio.cpp
        imageio.h
        linearSolver.cpp
        linearSolver.h
        Particle.cpp
        Particle.h
        RodConstraint.cpp
        RodConstraint.h
        Solver.cpp
        Solver.h
        SpringForce.cpp
        SpringForce.h
        TinkerToy.cpp
        MouseSpringForce.cpp
        MouseSpringForce.h
        AngularSpring.cpp
        AngularSpring.h
        JWJTranspose.cpp
        JWJTranspose.h
        Constraint.h)

add_executable(tue_scg_particle_systems ${SOURCE_FILES})

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(PNG REQUIRED)
include_directories( include ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} )

target_link_libraries(tue_scg_particle_systems ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${PNG_LIBRARY})