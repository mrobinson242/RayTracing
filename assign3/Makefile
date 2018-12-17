# Makefile 
# CSCI 420
# Assignment 3

# we assume the pic directory locates one level above,
# change PIC_PATH if this is not the case
PIC_PATH = $(abspath $(CURDIR)/../pic)
HEADER_PATH = $(abspath $(CURDIR)/headers)

INCLUDE = -I$(PIC_PATH)
HEADER = -I$(HEADER_PATH)
LIBRARIES = -L$(PIC_PATH) -framework OpenGL -framework GLUT -lpicio $(PIC_PATH)/libjpeg.a


COMPILER = g++
COMPILERFLAGS = -O3 $(INCLUDE)

PROGRAM = RayTracer
SOURCE = RayTracer.cpp
OBJECT = RayTracer.o VectorMath.o

.cpp.o: 
	$(COMPILER) -c $(COMPILERFLAGS) $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECT)
	$(COMPILER) $(COMPILERFLAGS) -o $(PROGRAM) $(OBJECT) $(LIBRARIES)

clean:
	-rm -rf core *.o *~ "#"*"#" $(PROGRAM)