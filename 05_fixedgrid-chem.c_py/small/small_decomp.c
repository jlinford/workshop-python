/*------------------------- BEGIN small_decomp.c BEGIN ------------------------*/
/* @file small_decomp.c                                                        */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.282053                                            */
/* @brief LU decomposition of the row-compressed sparse Jacobian               */
/*                                                                             */
/* LU decomposition of the row-compressed sparse Jacobian                      */
/*                                                                             */
/* This file was generated by Kppa: http://www.paratools.com/Kppa              */
/*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "small_parameters.h"
#include "small_sparse.h"
#include "small_decomp.h"


/*----------------------------------- Decomp ----------------------------------*/
/* In-place sparse LU decomposition                                            */
/*                                                                             */
/* @param[in,out] A Row-compressed matrix with zero fill                       */
/*-----------------------------------------------------------------------------*/
int Decomp(float A[JAC_LU_NZ])
{
  float t0;
  float t1;
  float t2;
  float t3;

  t0 = A[2]/A[0];
  A[2] = t0;
  A[4] = -A[1]*t0 + A[4];
  t1 = A[6]/A[0];
  A[6] = t1;
  A[8] = -A[1]*t1 + A[8];
  t2 = A[7]/A[3];
  A[7] = t2;
  A[8] = -A[4]*t2 + A[8];
  A[10] = A[10] - A[5]*t2;
  t3 = A[11]/A[3];
  A[11] = t3;
  A[12] = A[12] - A[4]*t3;
  A[14] = A[14] - A[5]*t3;
  t0 = A[12]/A[8];
  A[12] = t0;
  A[13] = A[13] - A[9]*t0;
  A[14] = -A[10]*t0 + A[14];
  t1 = A[15]/A[3];
  A[15] = t1;
  A[16] = A[16] - A[4]*t1;
  A[18] = A[18] - A[5]*t1;
  t2 = A[16]/A[8];
  A[16] = t2;
  A[17] = A[17] - A[9]*t2;
  A[18] = -A[10]*t2 + A[18];
  t3 = A[17]/A[13];
  A[17] = t3;
  A[18] = -A[14]*t3 + A[18];

  return 0;
}/* END Decomp */



/*--------------------------- END small_decomp.h END --------------------------*/