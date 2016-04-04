#****************************************************************************
#*			TAU Portable Profiling Package			   **
#*			http://www.cs.uoregon.edu/research/tau	           **
#****************************************************************************
#*    Copyright 1997  						   	   **
#*    Department of Computer and Information Science, University of Oregon **
#*    Advanced Computing Laboratory, Los Alamos National Laboratory        **
#****************************************************************************
#######################################################################
##                  pC++/Sage++  Copyright (C) 1993,1995             ##
##  Indiana University  University of Oregon  University of Rennes   ##
#######################################################################
 
CXX		= xlC_r 
CXX		= tau_cxx.sh

CC		= xlc_r 
CC		= tau_cc.sh

CFLAGS          = -O3  -openmp

LIBS            = 

LDFLAGS	      	= $(CFLAGS)

MAKEFILE      	= Makefile

PRINT	      	= pr

RM 	      	= /bin/rm -f

TARGET	      	= mandel

EXTRAOBJS     	= 

##############################################

mandel: 	mandel.o mytimer.o ppmwrite.o
	$(CXX) $(CFLAGS) mandel.o mytimer.o ppmwrite.o -o $@ $(LIBS)
mandel.o:	mandel.cpp
	$(CXX) $(CFLAGS) -c mandel.cpp
mytimer.o:	mytimer.c
	$(CC) $(CFLAGS) -c mytimer.c
ppmwrite.o:	ppmwrite.cpp
	$(CXX) $(CFLAGS) -c ppmwrite.cpp

clean: 	
	$(RM) mandel *.o  profile.* *.trc  *.ppm *.pomp* *.chk.* *opari.* *.pdb
##############################################
