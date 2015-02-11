PROJECT = TorusKnot
CXX = g++

LIB = libOGL.a

#HFILES = OGL_Vector.h OGL_Utility.h OGL_Camera.h local_classes.h global_prototypes.h
OFILES = OGL_Matrix.o OGL_Vector.o OGL_Utility.o OGL_Camera.o
AFILES = $(OFILES)

INCLUDES := -I ./ \

LINKS := \
	-L/usr/local/lib \

#CFLAGS = -g
UNAME	:= $(shell uname)
CFLAGS = -g -lm -lGL -lglut -lGLU -Wall -std=c++0x
OFLAGS = -c $(INCLUDES)

ifeq ($(UNAME), Darwin)
CFLAGS = -g -lm -framework OpenGL -framework GLUT
endif

all: $(OFILES) 
	ar rv $(LIB) $?

.C.o:
	${CXX} $(CFLAGS) $< $(OFLAGS)

.cpp:
	$(CXX) $(CFLAGS) $@.cpp $(INCLUDES) $(LIB) $(LINKS) -o $@

debug:
	make 'DFLAGS = /usr/lib/debug/malloc.o'

clean:
	rm *.o *~ ${PROJECT}

