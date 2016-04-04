

F90=tau_f90.sh -tau_options='-optVerbose -optMemDbg -optKeepFiles'

simple: simple.o
	@$(F90) simple.o -o simple
simple.o: simple.f90
	@$(F90) -c simple.f90
clean:
	/bin/rm -rf simple simple.o *.inst.* *.pdb  profile.* *.trc MULTI_* *.edf
	
	

