CC = tau_cc.sh
# Need to specify TAU_MAKEFILE environment variable. Or else uncomment the 
# following line on the laptop.
#CC = tau_cc.sh -tau_makefile=/usr/cta/pet/tau/i386_linux/lib/Makefile.tau-icpc-pthread-pdt

threaded: threaded.o
	$(CC) -o threaded threaded.o -lpthread
threaded.o: threaded.c
	$(CC) -c threaded.c
clean:
	/bin/rm -f threaded.o threaded profile.* *.elg *.trc *.edf *.pdb *.inst.c

