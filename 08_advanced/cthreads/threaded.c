#include <stdio.h>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>

int fourth(void)
{
  printf("Reached fourth\n");

  sleep(3);
  return 0;
}

int third(void)
{
  printf("third calling fourth()\n");
  fourth();
  return 0;
}

int second(void)
{
  printf("second calling third\n");
  third();
  return 0;
}

int first(void)
{ 
  printf("first.. calling second \n");
  second();
  return 0;
}

int work (void)
{
  printf("Hello. This is thread work calling first\n");
  sleep(5);
  first(); 
  return 0;
}
void * threaded_func(void *data)
{
  work(); /* work done by this thread */
  return NULL;
}

int main (int argc, char **argv)
{
  int ret, i;
  pthread_attr_t  attr;
  pthread_t	  tid;

  pthread_attr_init(&attr);
  
  printf("Started Main...\n");

  if (ret = pthread_create(&tid, NULL, threaded_func, NULL) ) 
  { 
    printf(" pthread_create fails ret = %d\n", ret); 
    exit(1);
  }


  if (ret = pthread_join(tid, NULL) ) 
  {
    printf(" pthread_join error  ret = %d\n", ret);
    exit(1); 
  }

  /* prior to exiting, print statistics related to user defined events */
  printf("Exiting main...\n");
  return 0;
}
/***************************************************************************
 * $RCSfile: hello.c,v $   $Author: sameer $
 * $Revision: 1.6 $   $Date: 1999/06/20 05:01:24 $
 * POOMA_VERSION_ID: $Id: hello.c,v 1.6 1999/06/20 05:01:24 sameer Exp $
 ***************************************************************************/

