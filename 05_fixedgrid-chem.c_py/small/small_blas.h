/*-------------------------- BEGIN small_blas.h BEGIN -------------------------*/
/* @file small_blas.h                                                          */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.220711                                            */
/* @brief Basic linear algebra subprogram definitions                          */
/*                                                                             */
/* A reduced set of BLAS routines optimized for Kppa-generated solvers         */
/*                                                                             */
/* This file was generated by Kppa: http://www.paratools.com/Kppa              */
/*-----------------------------------------------------------------------------*/

#ifndef __SMALL_BLAS_H__
#define __SMALL_BLAS_H__



#ifdef __cplusplus
extern "C" {
#endif


/*----------------------------------- WCOPY -----------------------------------*/
/* Copies vector x to vector y: y <= x                                         */
/* Like the BLAS {S,D}COPY(N,X,1,Y,1)                                          */
/*                                                                             */
/* @param[in]     n Vector length                                              */
/* @param[in]     x Vector x                                                   */
/* @param[out]    y Vector y                                                   */
/*-----------------------------------------------------------------------------*/
void WCOPY(size_t const n, float const  x[/* n */], float y[/* n */]);


/*----------------------------------- WSCAL -----------------------------------*/
/* Constant times a vector: x <= alpha*x                                       */
/* Like the BLAS {S,D}SCAL(N,alpha,X,1)                                        */
/*                                                                             */
/* @param[in]     n     Vector length                                          */
/* @param[in]     alpha Scalar                                                 */
/* @param[in,out] x     Vector x                                               */
/*-----------------------------------------------------------------------------*/
void WSCAL(size_t const n, float const alpha, float x[/* n */]);


/*----------------------------------- WAXPY -----------------------------------*/
/* Constant times a vector plus a vector: y <= y + alpha*x                     */
/* Like the BLAS {S,D}AXPY(N,alpha,X,1,Y,1)                                    */
/*                                                                             */
/* @param[in]     n     Vector length                                          */
/* @param[in]     alpha Scalar                                                 */
/* @param[in]     x     Vector x                                               */
/* @param[in,out] y     Vector y                                               */
/*-----------------------------------------------------------------------------*/
void WAXPY(size_t const n, float const alpha, float const  x[/* n */], float
    y[/* n */]);


/*----------------------------------- WYMXDA ----------------------------------*/
/* Difference of two vectors divided by a constant: z <= (y - x) / alpha       */
/*                                                                             */
/* @param[in]     n     Vector length                                          */
/* @param[in]     x     Vector x                                               */
/* @param[in]     y     Vector y                                               */
/* @param[in]     alpha Scalar                                                 */
/* @param[out]    z     Vector z                                               */
/*-----------------------------------------------------------------------------*/
void WYMXDA(size_t const n, float const  x[/* n */], float const  y[/* n */],
    float const alpha, float z[/* n */]);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SMALL_BLAS_H__ */
/*---------------------------- END small_blas.h END ---------------------------*/
