/*----------------------- BEGIN small_rosenbrock.h BEGIN ----------------------*/
/* @file small_rosenbrock.h                                                    */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.349747                                            */
/* @brief Solves the system y' = F(t,y) using a Rosenbrock method              */
/*                                                                             */
/* Solves the system y' = F(t,y) using a Rosenbrock method defined by:         */
/*                                                                             */
/*     G = 1 / (H*gamma) - Jacobian(t0,Y0)                                     */
/*     T_i = t0 + Alpha(i) * H                                                 */
/*     Y_i = Y0 + \sum_{j=1}^{i-1} A(i,j)*K_j                                  */
/*     G * K_i = F(T_i, Y_i) + \sum_{j=1}^S C(i,j)/H * K_j                     */
/*               + gamma(i)*dF/dT(t0, Y0)                                      */
/*     Y1 = Y0 + \sum_{j=1}^S M(j)*K_j                                         */
/*                                                                             */
/* For details on Rosenbrock methods and their implementations:                */
/*     (1) E. Harier and G. Wanner,                                            */
/*         "Solving Ordenary Differential Equations II: stiff and              */
/*         differential-algebraic problems." Computational Mathematics,        */
/*         Springer-Verlag (1996)                                              */
/*     (2) KPP - the Kinetic PreProcessor.                                     */
/*         http://people.cs.vt.edu/~asandu/Software/Kpp/                       */
/*                                                                             */
/* Rosenbrock implementations in both (1) and (2) inspired this code.          */
/* This code presents an interface similar to the KPP implementation           */
/* for compatibility with existing systems.                                    */
/*                                                                             */
/* -- Explanation of integer input parameters:                                 */
/*                                                                             */
/*     idata[0] == 0 : F = F(t,y) Depends on T (non-autonomous).               */
/*              != 0 : F = F(y)   Independent of T (autonomous).               */
/*     idata[1] == 0 : Use all values in tolerance vectors.                    */
/*              != 0 : Use only the first value in the tolerance vectors.      */
/*     idata[2] == 0 : Maximum number of integration steps = 100000.           */
/*              != 0 : Maximum number of integration steps = idata[2].         */
/*     idata[3] == 0 : Method is Ros4.                                         */
/*              == 1 : Method is Ros2.                                         */
/*              == 2 : Method is Ros3.                                         */
/*              == 3 : Method is Ros4.                                         */
/*              == 4 : Method is Rodas3.                                       */
/*              == 5 : Method is Rodas4.                                       */
/*              >= 6 : Error.                                                  */
/*     idata[4] == 0 : Assume tolerance vectors are reasonably valued.         */
/*              != 0 : Check tolerance vectors for unreasonable values.        */
/*                                                                             */
/* -- Explanation of real value input parameters:                              */
/*                                                                             */
/*     rdata[0]: Lower bound on the integration step size.                     */
/*               Default: 0.0                                                  */
/*     rdata[1]: Upper bound on the integration step size.                     */
/*               Default: abs(tend - tstart)                                   */
/*     rdata[2]: Starting value for the integration step size.                 */
/*               Default: minimum step size                                    */
/*     rdata[3]: Lower bound on step decrease factor.                          */
/*               Default: 0.2                                                  */
/*     rdata[4]: Upper bound on step increase factor.                          */
/*               Default: 6.0                                                  */
/*     rdata[5]: Step decrease factor after step rejection.                    */
/*               Default: 0.1                                                  */
/*     rdata[6]: Safety factor in computation of new step size.                */
/*               Default: 0.9                                                  */
/*                                                                             */
/* -- Explanation of integer output parameters:                                */
/*                                                                             */
/*     idata[10]: Number of function evaluations.                              */
/*     idata[11]: Number of Jacobian evaluations.                              */
/*     idata[12]: Number of steps taken.                                       */
/*     idata[13]: Number of accepted steps.                                    */
/*     idata[14]: Number of rejected steps.                                    */
/*     idata[15]: Number of LU decompositions.                                 */
/*     idata[16]: Number of forward/backward substitutions.                    */
/*     idata[17]: Number of singular matrix decompositions.                    */
/*     idata[19]: Integrator exit status.                                      */
/*                Zero indicates success.                                      */
/*                Positive values indicate success with warning.               */
/*                Negative values indicate failure.                            */
/*                                                                             */
/* -- Explanation of real-value output parameters:                             */
/*                                                                             */
/*     rdata[10]: The time corresponding to the computed Y upon return.        */
/*     rdata[11]: The last accepted step before exit.                          */
/*                Use this value as rdata[2] in subsequent runs.               */
/*     rdata[12]: Scaled norm of the error vector on exit.                     */
/*                                                                             */
/* This file was generated by Kppa: http://www.paratools.com/Kppa              */
/*-----------------------------------------------------------------------------*/

#ifndef __SMALL_ROSENBROCK_H__
#define __SMALL_ROSENBROCK_H__



#ifdef __cplusplus
extern "C" {
#endif


/*--------------------------------- Integrate ---------------------------------*/
/* Kppa-generated time stepping integrator                                     */
/*                                                                             */
/* @param[in,out] var    Variable species concentrations                       */
/* @param[in,out] fix    Fixed species concentrations                          */
/* @param[in]     idx    Current grid cell index                               */
/* @param[in]     tstart Integration start time                                */
/* @param[in]     tend   Integration end time                                  */
/* @param[in]     abstol Absolute integration tolerances for variable species  */
/* @param[in]     reltol Relative integration tolerances for variable species  */
/* @param[in,out] idata  Integer integration in/out parameters                 */
/* @param[in,out] rdata  Real value integration in/out parameters              */
/*-----------------------------------------------------------------------------*/
void Integrate(float var[5], float fix[2], size_t const idx, float const
    tstart, float const tend, float const  abstol[5], float const  reltol[5],
    int idata[20], float rdata[20]);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* __SMALL_ROSENBROCK_H__ */
/*------------------------- END small_rosenbrock.h END ------------------------*/
