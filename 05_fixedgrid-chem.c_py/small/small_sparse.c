/*------------------------- BEGIN small_sparse.c BEGIN ------------------------*/
/* @file small_sparse.c                                                        */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.168686                                            */
/* @brief Data and utilities for row-compressed sparse matrices                */
/*                                                                             */
/* The following matrices are represented here in row-compressed form:         */
/* @li The Jacobian                                                            */
/* @li The LU decomposition of the Jacobian                                    */
/* @li The stoichiometric matrix                                               */
/*                                                                             */
/* This file was generated by Kppa: http://www.paratools.com/Kppa              */
/*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "small_sparse.h"


/* Number of nonzero entries in the Jacobian */
#define JAC_NZ           ((size_t)(19))

/* Number of nonzero entries in the LU decomposition of the Jacobian */
#define JAC_LU_NZ        ((size_t)(19))

/* Number of nonzero entries in the stoichiometric matrix */
#define STOICH_NZ        ((size_t)(34))

/* Number of nonzero entries in the left-side stoichiometric matrix */
#define LHS_STOICH_NZ    ((size_t)(20))

/* Number of nonzero entries in the right-side stoichiometric matrix */
#define RHS_STOICH_NZ    ((size_t)(16))


/* Row indices of elements in the row-compressed Jacobian */
int const JAC_IROW[19] = { 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4,
    4, 4 };

/* Column indices of elements in the row-compressed Jacobian */
int const JAC_ICOL[19] = { 0, 2, 0, 1, 2, 4, 0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2,
    3, 4 };

/* Start-of-row indices in the row-compressed Jacobian */
int const JAC_CROW[6] = { 0, 2, 6, 11, 15, 19 };

/* Diagonal indices in the row-compressed Jacobian */
int const JAC_DIAG[6] = { 0, 3, 8, 13, 18, 19 };

/* Row indices of elements in the row-compressed LU decomposition of the Jacobian */
int const JAC_LU_IROW[19] = { 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4,
    4, 4, 4 };

/* Column indices of elements in the row-compressed LU decomposition of the Jacobian */
int const JAC_LU_ICOL[19] = { 0, 2, 0, 1, 2, 4, 0, 1, 2, 3, 4, 1, 2, 3, 4, 1,
    2, 3, 4 };

/* Start-of-row indices in the row-compressed LU decomposition of the Jacobian */
int const JAC_LU_CROW[6] = { 0, 2, 6, 11, 15, 19 };

/* Diagonal indices in the row-compressed LU decomposition of the Jacobian */
int const JAC_LU_DIAG[6] = { 0, 3, 8, 13, 18, 19 };

/* Row indices of elements in the row-compressed stoichiometric matrix */
int const STOICH_IROW[34] = { 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    3, 3, 3, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7 };

/* Column indices of elements in the row-compressed stoichiometric matrix */
int const STOICH_ICOL[34] = { 4, 5, 6, 0, 1, 2, 3, 5, 8, 9, 1, 2, 3, 4, 6, 7,
    7, 8, 9, 7, 8, 9, 0, 1, 2, 3, 4, 6, 7, 8, 0, 2, 4, 9 };

/* Start-of-row indices in the row-compressed stoichiometric matrix */
int const STOICH_CROW[9] = { 0, 3, 10, 16, 19, 22, 22, 30, 34 };

/* Stoichiometric coefficients */
float const STOICH[34] = { 1.0, -1.0, -1.0, 2.0, -1.0, 1.0, -1.0, 1.0, -1.0,
    1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0,
    -1.0, -1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0 };

/* Row indices of elements in the row-compressed left-side stoichiometric matrix */
int const LHS_STOICH_IROW[20] = { 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 4, 4, 5,
    6, 6, 7, 7, 7, 7 };

/* Column indices of elements in the row-compressed left-side stoichiometric matrix */
int const LHS_STOICH_ICOL[20] = { 5, 6, 1, 3, 8, 2, 3, 4, 6, 7, 7, 8, 9, 5,
    0, 1, 0, 2, 4, 9 };

/* Start-of-row indices in the row-compressed left-side stoichiometric matrix */
int const LHS_STOICH_CROW[9] = { 0, 2, 5, 10, 11, 13, 14, 16, 20 };

/* Left-side stoichiometric coefficients */
float const LHS_STOICH[20] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

/* Row indices of elements in the row-compressed right-side stoichiometric matrix */
int const RHS_STOICH_IROW[16] = { 0, 1, 1, 1, 1, 2, 3, 3, 4, 5, 6, 6, 6, 6,
    6, 6 };

/* Column indices of elements in the row-compressed right-side stoichiometric matrix */
int const RHS_STOICH_ICOL[16] = { 4, 0, 2, 5, 9, 1, 8, 9, 7, 5, 2, 3, 4, 6,
    7, 8 };

/* Start-of-row indices in the row-compressed right-side stoichiometric matrix */
int const RHS_STOICH_CROW[9] = { 0, 1, 5, 6, 8, 9, 10, 16, 16 };

/* Right-side stoichiometric coefficients */
float const RHS_STOICH[16] = { 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0 };


/*------------------------------------ CSR ------------------------------------*/
/* Retrieves an element from a compressed sparse row matrix                    */
/*                                                                             */
/* @param[in]     nz   Matrix nonzero values                                   */
/* @param[in]     crow Row start indices                                       */
/* @param[in]     icol Column indices                                          */
/* @param[in]     row  Row index into matrix                                   */
/* @param[in]     col  Column index into matrix                                */
/*-----------------------------------------------------------------------------*/
float CSR(float const * nz, int const * crow, int const * icol, int const
    row, int const col)
{
  /* i-loop index */
  int i;

  for(i = crow[row]; i < crow[row + 1]; ++i) {
    if (icol[i] == col) {
      return nz[i];
    } /* END for(icol[i] == col) */
  } /* END for i=crow[row]:crow[row + 1] */
  return 0;
}/* END CSR */



/*--------------------------- END small_sparse.h END --------------------------*/
