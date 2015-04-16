/*----------------------- BEGIN small_integrate.c BEGIN -----------------------*/
/* @file small_integrate.c                                                     */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.356016                                            */
/* @brief Interface to time stepping integrator                                */
/*                                                                             */
/* Definitions of interface functions for the Kppa-generated                   */
/* time stepping integrator.  These are the Kppa "entry point" routines.       */
/*                                                                             */
/* This file was generated by Kppa: http://www.paratools.com/Kppa              */
/*-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "small_parameters.h"
#include "small_rosenbrock.h"
#include "small_integrate.h"


/*------------------------------- GridIntegrate -------------------------------*/
/* Applies the Kppa-generated integrator to the grid                           */
/*                                                                             */
/* @param[in]     ncells Number of grid cells                                  */
/* @param[in,out] conc   Species concentrations                                */
/* @param[in]     tstart Integration start time                                */
/* @param[in]     tend   Integration end time                                  */
/* @param[in]     abstol Absolute integration tolerances for variable species  */
/* @param[in]     reltol Relative integration tolerances for variable species  */
/* @param[in,out] idata  Integer integration in/out parameters                 */
/* @param[in,out] rdata  Real value integration in/out parameters              */
/* @param[out]    lastH  Last timestep in each grid cell                       */
/*-----------------------------------------------------------------------------*/
int GridIntegrate(size_t const ncells, float conc[], float const tstart,
    float const tend, float const  abstol[NVAR], float const  reltol[NVAR],
    int idata[20], float rdata[20], float lastH[/* ncells */])
{
    /* Return value */
    int retval = 0;
    /* Iterators */
    int i, j;
    /* Pointer to variable concentrations in grid cell */
    float * var;
    /* Pointer to fixed concentrations in grid cell */
    float * fix;

    #pragma omp parallel for \
                default(shared) \
                private(i,j,var,fix) \
                firstprivate(idata, rdata) \
                lastprivate(idata, rdata) \
                reduction(|:retval)
    for(i=0; i<ncells; ++i) {
        /* Get concentrations for the current grid cell */
        var = conc + i*NSPEC;
        fix = var + NVAR;

        /* Invoke the integrator */
        Integrate(var, fix, i, tstart, tend, abstol, reltol, idata, rdata);

        /* Save the last timestep for future use */
        if(lastH) {
            lastH[i] = rdata[11];
        }

        /* Process integrator return code */
        if (idata[19] < 0) {
            printf("Kppa: CELL %d -- INTEGRATION FAILED\n", i);
            #pragma omp critical
            {
                for(j=0; j<20; ++j)
                    printf("Kppa: CELL %d -- idata[%d] = %d\n", i, j, idata[j]);
                for(j=0; j<20; ++j)
                    printf("Kppa: CELL %d -- rdata[%d] = %g\n", i, j, rdata[j]);
            }
            if (idata[19] < retval)
                retval = idata[19];
        } else if (idata[19] > 0) {
            printf("Kppa: CELL %d -- INTEGRATION COMPLETED WITH WARNING\n", i);
            if (retval >= 0 && idata[19] > retval)
                retval = idata[19];
        }
    }
    
    return retval;
}/* END GridIntegrate */


/*------------------------- END small_integrate.h END -------------------------*/
