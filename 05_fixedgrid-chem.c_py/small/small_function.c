/*------------------------ BEGIN small_function.c BEGIN -----------------------*/
/* @file small_function.c                                                      */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.255338                                            */
/* @brief The ODE function of the chemical model                               */
/*                                                                             */
/* The ODE function of the chemical model                                      */
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
#include "small_function.h"


/*------------------------------------ Fun ------------------------------------*/
/* The ODE function of the chemical model                                      */
/*                                                                             */
/* @param[in]     var    Variable species concentrations                       */
/* @param[in]     fix    Fixed species concentrations                          */
/* @param[in]     rct    Reaction rates                                        */
/* @param[out]    vardot The ODE function                                      */
/*-----------------------------------------------------------------------------*/
void Fun(float const  var[NVAR], float const  fix[NFIX], float const
    rct[NREACT], float vardot[NVAR])
{
  float r0;
  float r1;
  float r2;
  float r3;
  float r4;
  float r5;
  float r6;
  float r7;
  float r8;
  float r9;

  r0 = fix[1]*rct[0];
  r1 = 8.018e-17*fix[1]*var[1];
  r2 = rct[2]*var[2];
  r3 = 1.576e-15*var[1]*var[2];
  r4 = rct[4]*var[2];
  r5 = 7.11e-11*fix[0]*var[0];
  r6 = 1.2e-10*var[0]*var[2];
  r7 = 6.062e-15*var[2]*var[3];
  r8 = 1.069e-11*var[1]*var[4];
  r9 = rct[9]*var[4];

  vardot[0] = r4 - r5 - r6;
  vardot[1] = 2.0*r0 - r1 + r2 - r3 + r5 - r8 + r9;
  vardot[2] = r1 - r2 - r3 - r4 - r6 - r7;
  vardot[3] = -r7 + r8 + r9;
  vardot[4] = r7 - r8 - r9;

}/* END Fun */



/*-------------------------- END small_function.h END -------------------------*/
