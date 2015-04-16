"""
TAU wrapper for fixedgrid

To run with TAU:
mpirun -np 4 tau_exec -T python,mpi python wrapper.py 
"""
import tau
tau.run('import fixedgrid')

