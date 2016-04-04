F90=mpif90
#F90=tau_f90.sh
FFLAGS= -O3 -g
LIBS=

matmult: matmult.f90
	$(F90) $(FFLAGS) matmult.f90 -o matmult $(LIBS)
clean:
	/bin/rm -rf *.o matmult profile.* MULTI__* *.trc *.edf *.otf *.def* *events*  *.out *.err
