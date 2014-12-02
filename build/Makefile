# A Generic Makefile for ADMB programs that also includes additonal libraries.
# Developed for Mac OSx using the clang++ compiler
# Author: Steve Martell & John Sibert

# establish the C++ compiler (on Mac OSX use clang++)
CC=gcc
CXX=g++
#CC=clang
#CXX=clang++
# and linker
LL = $(CC)
LD = $(CXX)
# Remove macro
RM=rm -fv

# identify some extra file name suffixes
#.SUFFIXES: .tpl .cpp .o 
.SUFFIXES:

# tell make not to delete these intermediate targets
.PRECIOUS: %.c %.cpp %.o

# make some special PHONY targets
.PHONY: all help rules clean dox



ifndef ADMB_HOME
  ADMB_HOME=/usr/local/admb
endif

ADMB_XML=/home/jsibert/admb/trunk/tests/xml
TPL_DIR = ../src
SRC_DIR = ../src
VPATH=$(TPL_DIR):$(SRC_DIR):$(ADMB_XML)


# set up ADMB flags and appropriate libraries
# make the "safe" version by default
# to make "the optimized" version, type  `make OPT=TRUE
ifeq ($(OPT),TRUE)
  CC_OPT = -O3 -DOPT_LIB
  LFLAGS= $(ADMB_HOME)/lib/libadmbo.a $(ADMB_HOME)/lib/libadmb-contribo.a -lxml2
else
  CC_OPT = -ggdb -O0
  LFLAGS= $(ADMB_HOME)/lib/libadmb.a $(ADMB_HOME)/lib/libadmb-contrib.a -lxml2
  OPT = FALSE
endif

CFLAGS = -m64 -Wall $(CC_OPT) -I../tpl -I../src -I$(ADMB_HOME)/include -I$(ADMB_HOME)/include/contrib -I$(ADMB_HOME)/contrib/gdbprintlib -I$(ADMB_XML) -I/usr/include/libxml2 -D__GNUDOS__ -Dlinux -DUSE_LAPLACE
 
CXXFLAGS=$(CFLAGS)
 
BIN_DIR = .
EXE =
LL=$(CC)

OBJECTS=ADMB_XMLDoc.o model_xml.o

default: xssams
#all:\
#read_plot_rep \
#yft0


:%.o
	@echo ------- Link $^; 
	$(LD) -o ../run/$@ $^ $(LFLAGS)

%.o:%.cpp
	@echo ------- Compile $<; 
	$(CXX) -c $(CFLAGS) -o $@ $<

%.cpp:%.tpl 
	@echo ------- Translate $<; 
	$(ADMB_HOME)/bin/tpl2rem $(basename $<)

xssams: xssams.o
	@echo ------- Link $^; 
	$(LD) -o ../run/$@ $^ $(LFLAGS)

#xssams.cpp: xssams.tpl
#	@echo ------- Translate $<; 
#	$(ADMB_HOME)/bin/tpl2rem $(basename $<)
#
#libxfer.a : $(OBJECTS)
#	@echo ------- Create $@; 
#	@echo         from $(OBJECTS)
#	ar -vrc $@ $^
#
#read_plot_rep: read_plot_rep.o libxfer.a
#	@echo ------- Link $@; 
#	$(LD) -o $@ $^ $(LFLAGS) -L. -lxfer
#
#yft0: yft0.o libxfer.a
#	@echo ------- Link $^; 
#	$(LD) -o $@ $^ $(LFLAGS) -L. -lxfer

#####################################################################

REV = `svnversion .`
ADMB_VERSION=`cat $(ADMB_HOME)/VERSION`
rules:
	@echo VPATH = $(VPATH)
	@echo vpath = $(vpath)
	@echo OPT = $(OPT)
	@echo CC_OPT = $(CC_OPT)
	@echo CURDIR = $(CURDIR)
	@echo PWD = $(PWD)
	@echo SRC_DIR=$(SRC_DIRS)
	@echo OBJ_DIR = $(OBJ_DIR)
	@echo OSTYPE = $(OSTYPE)
	@echo OS     = $(OS)
	@echo TERM = $(TERM)
	@echo BIN_DIR = $(BIN_DIR)
	@echo ADMB_HOME = $(ADMB_HOME)
	@echo ADMB_VERSION = $(ADMB_VERSION)
	@echo CC = $(CC)
	$(CC) --version
	@echo CXX = $(CXX)
	$(CXX) --version
	@echo LL = $(LL)
	@echo CFLAGS = $(CFLAGS)
	@echo CXXFLAGS = $(CXXFLAGS)
	@echo LFLAGS = $(LFLAGS)
	@echo EXE = $(EXE)

clean:
	@rm -vf ./*.o
	@rm -vf $(TPL_DIR)/*.cpp
	@rm -vf $(TPL_DIR)/*.htp
	@rm -vf libxfer.a
	@rm -vf yft0
	@rm -vf read_plot_rep
	



