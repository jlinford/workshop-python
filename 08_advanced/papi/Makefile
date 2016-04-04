CC=tau_cc.sh -tau_options='-optTauSelectFile=select.tau -optVerbose -optKeepFiles '

#CFLAGS = -O3
APP=loop_test
$(APP): $(APP).o
	$(CC) $(APP).o -o $@ 
$(APP).o: $(APP).c
	$(CC) -c $(APP).c  $(CFLAGS)
# -UTAU_MPI just protects the compilation in case you use an MPI enabled 
# stub makefile with a sequential program
clean:
	/bin/rm -rf *.o *.pdb *.inst.* $(APP) profile.* MULTI* run.job.*
