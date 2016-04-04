#F90 = mpif90
F90 = tau_f90.sh

LIBS = 

ring: ring.o
	$(F90) ring.o -o ring $(LIBS)

ring.o: ring.f90
	$(F90) -c ring.f90

clean:
	/bin/rm -f ring.o ring *.trc profile.* *.edf *.pdb *.inst.* *.slog2 *.otf* *.def *.events
