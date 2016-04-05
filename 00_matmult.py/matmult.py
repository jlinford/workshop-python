"""
Matrix multiplication examples
"""
import random

MATRIX_SIZE = 200
NRA = MATRIX_SIZE
NCA = MATRIX_SIZE
NCB = MATRIX_SIZE

def naive_matrix(m, n):
  """
  Stupidly creates an MxN matrix of random data using Python lists
  """
  return [random.sample(xrange(n), n) for _ in xrange(m)]

def naive_matmult(c, a, b):
  """
  Stupidly multiply list-based matrices
  """
  rows_a = len(a)
  cols_a = len(a[0])
  cols_b = len(b[0])
  for i in xrange(rows_a):
    for j in xrange(cols_b):
      for k in xrange(cols_a):
        c[i][j] += a[i][k] * b[k][j]

def naive_example():
  """
  Multiply matrices in the worst way possible
  """
  print 'Multiplying matrices in a very stupid way'
  a = naive_matrix(NRA, NCA)
  b = naive_matrix(NCA, NCB)
  c = naive_matrix(NRA, NCB)
  naive_matmult(c, a, b)

def numpy_example():
  """
  Multiply matrices via numpy
  """
  print 'Multiplying matrices via numpy'
  import numpy as np
  a = np.random.rand(NRA, NCA)
  b = np.random.rand(NCA, NCB)
  c = np.random.rand(NRA, NCB)
  c = np.dot(a, b)
  
if __name__ == "__main__":
  naive_example()
  numpy_example()
