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
.SUFFIXES: .tpl

# tell make not to delete these intermediate targets
.PRECIOUS: %.c %.cpp %.o

# make some special PHONY targets
.PHONY: all clean
all: strat csm


strat: strat.cpp
	g++ -m64 -Wall -ggdb -DTEST_CODE -DSAFE_ARRAYS -Isrc -I${ADMB_HOME}/include -I${ADMB_HOME}/contrib/include -I${JNIGRAPHICS_PATH} -I${JAVA_HOME}/include -I${JAVA_HOME}/include/${OS} -D__GNUDOS__ -Dlinux -DUSE_LAPLACE -UOLD_DO_KLUDGE -c strat.cpp -o strat.o
	g++ strat.o ${ADMB_HOME}/lib/libadmb.a ${ADMB_HOME}/lib/libadmb-contrib.a -o strat
 
csm: csm.tpl
	admb -g csm


clean:
	@rm -vf *.o
	@rm -vf *.obj
	@rm -vf *.p??
	@rm -vf *.r??
	@rm -vf *.b??
	@rm -vf *.eva
	@rm -vf *.cor
	@rm -vf *.std
	@rm -vf *.htp
	@rm -vf *.log
	@rm -vf *.def
	@rm -vf admodel.*
	@rm -vf variance
	@rm -vf strat
	@rm -vf csm
	@rm -vf csm.cpp

	



