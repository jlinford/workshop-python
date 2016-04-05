#!/bin/tcsh

setenv TAU_TRACK_SIGNALS 1

#run the executable
mpirun -np 4 tau_python ./samarcrun.py

