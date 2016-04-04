#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>

int foo(int x) {
  sleep(x);
  return 0;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  TAU_START("sleeping for 4 secs");
  foo(4);
  TAU_STOP("sleeping for 4 secs");

  MPI_Finalize();
}

