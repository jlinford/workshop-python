#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <stdlib.h>

struct node {
  int id;
  struct node *next;
};

int bar(int x) {
  int y;
  struct node *t = (struct node *)malloc(sizeof(struct node));
  t->next = NULL;
  t->id = x;
  printf("t -> id = %d\n", t->id);
  y = t->next->id;
  printf("y=%d\n",y);
  return x;
}

int foo(int x) {
  printf("foo: x = %d\n", x);
  bar(x);
  return x;
}

int main(int argc, char **argv) {
  int ret;
  MPI_Init(&argc, &argv);
  ret = foo(29);
  MPI_Finalize();
  return ret;
}
