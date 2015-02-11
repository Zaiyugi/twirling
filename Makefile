
.SUFFIXES: .cpp .C .h .o

DATE := $(shell date +%F)
UNAME := $(shell uname)

ROOTDIR = .
LIB_DIR = ./lib
SRC_DIR = src
OBJ_DIR = build
INCLUDE = include
HT := /usr
HDSO := /usr/lib64

VPATH = src:include

OFILES = Matrix.o LinearAlgebra.o

OBJS = $(patsubst %, $(OBJ_DIR)/%, $(OFILES))

GL_DIR = $(ROOTDIR)/opengl

VR_LIB = $(LIB_DIR)/libVR.a

OPENMESH_DIR := /usr/local/OpenMesh/include
OPENMESH_LIB := /usr/local/OpenMesh/lib/OpenMesh
OPENMESH_CORE := OpenMeshCore

ATB_DIR := $(HOME)/local_lib/AntTweakBar/include
ATB_LIB := $(HOME)/local_lib/AntTweakBar/lib

LINKS := \
	-L$(LIB_DIR) \
	-L/usr/local/lib \
	-L$(HDSO) \
	-L$(ATB_DIR) \
	-L$(OPENMESH_LIB) -l$(OPENMESH_CORE) \
	-lVR \
	-lOGL \

INCLUDES := \
	-I. \
	-I./include \
	-I$(GL_DIR) \
	-I$(OPENMESH_DIR) \
	-I$(ATB_DIR) \

CXX = clang++
CFLAGS = -g -O2 -D_THREAD_SAFE -pthread -std=c++0x
OFLAGS = -c $(INCLUDES) 
GLFLAGS = -lAntTweakBar -lGL -lglut -lGLU

all: $(OBJS)
	ar rv $(VR_LIB) $?
	make gl

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C $(INCLUDE)/%.h
	$(CXX) $< -Wall $(CFLAGS) $(OFLAGS) -o $@

gl:
	cd $(GL_DIR); make all

launch: horn_modeling.cpp
	$(CXX) horn_modeling.cpp $(CFLAGS) $(GLFLAGS) $(INCLUDES) $(LINKS) -o launch

print:
	$@

clean:
	rm $(OBJ_DIR)/*.o $(LIB_DIR)/*.a launch

