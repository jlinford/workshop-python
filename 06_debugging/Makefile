#Makefile for driver directory of SAMARC

#### Defns ##########
CC            = mpicc
CXX           = mpicxx
F90           = mpif90
MPIF90        = mpif90

CFLAGS        = -g
CXXFLAGS      = -fPIC -fpermissive -g
FFLAGS        = -g 

PYTHON_PREFIX = /opt/ParaTools/anaconda2/4.0.0
NUMPY_PREFIX = $(PYTHON_PREFIX)/lib/python2.7/site-packages/numpy
PYINC = -I$(PYTHON_PREFIX)/include/python2.7 -I$(NUMPY_PREFIX)/core/include
PYLIB = -L$(PYTHON_PREFIX)/lib -lpython2.7

#### Codes ##########
CXX_OBJS  = SAMINT.o samint_wrap.o

F90_OBJS  = fortmthds.o

PY_OBJS   = pycintfc.o

#### Rules ##########
all:	$(CXX_OBJS) $(F90_OBJS) $(PY_OBJS)
	$(CXX) -fPIC -shared $(PY_OBJS) $(CXX_OBJS) $(F90_OBJS) \
	-o _samint.so $(PYLIB)

pycintfc.o: pycintfc.C
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c pycintfc.C $(PYINC) 

samint_wrap.o: samint_wrap.c
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c samint_wrap.c $(PYINC) 

samint_wrap.c: samint.i
	swig -python -I../  -I$(SWIG_LIB) -I$(SWIG_LIB)/python samint.i

fortmthds.o: fortmthds.f
	$(F90) $(FFLAGS) -c fortmthds.f


#### clean ##########
cleanc:
	$(RM) *.o 

clean:
	$(RM) _samint.so
	$(RM) *.o *.pyc *~

