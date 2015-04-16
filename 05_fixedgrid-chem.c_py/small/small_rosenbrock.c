/*----------------------- BEGIN small_rosenbrock.c BEGIN ----------------------*/
/* @file small_rosenbrock.c                                                    */
/* @author jlinford                                                            */
/* @date 2015-04-15 21:17:47.347415                                            */
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
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "small_parameters.h"
#include "small_sparse.h"
#include "small_blas.h"
#include "small_rates.h"
#include "small_function.h"
#include "small_decomp.h"
#include "small_solve.h"
#include "small_jacobian.h"
#include "small_rosenbrock.h"


#include <float.h>

/* Minimum time delta */
#define MIN_DELT ((float)10.0 * FLT_EPSILON)



/*----------------------------------------------------------------------------*/
/* A two-stage L-stable method of order 2                                     */
/*                                                                            */
/* E. Harier and G. Wanner, "Solving Ordenary Differential Equations II:      */
/* stiff and differential-algebraic problems." Computational Mathematics,     */
/* Springer-Verlag (1996)                                                     */
/*                                                                            */
/* @param[out] name     Method name                                           */
/* @param[out] nStage   Number of method stages                               */
/* @param[out] invLoEst One divided by the estimation of local order          */
/* @param[out] M        Coefficients for new step solution                    */
/* @param[out] E        Coefficients for error estimation                     */
/* @param[out] A        Lower triangular coefficient matrix                   */
/* @param[out] C        Lower triangular coefficient matrix                   */
/* @param[out] alpha    Y at stage i is approx. Y(T + H*Alpha_i)              */
/* @param[out] gamma    Stage i Gamma = sum(gamma[j])                         */
/* @param[out] F        Function evaluation flags                             */
/*----------------------------------------------------------------------------*/
void InitRos2(char ** name, int * nStage, float * invLoEst, float M[], float E[],
        float A[], float C[], float alpha[], float gamma[], char F[])
{
    *name = "Ros2";

    *nStage = 2;

    *invLoEst = 0.5; /* 1 / 2 */

    M[0] = 0.8786796564403575;
    M[1] = 0.2928932188134525;

    E[0] = 0.2928932188134525;
    E[1] = 0.2928932188134525;

    A[0] = 0.585786437626905;

    C[0] = -1.17157287525381;

    alpha[0] = 0.0;
    alpha[1] = 1.0;

    gamma[0] =  1.7071067811865475;
    gamma[1] = -1.7071067811865475;

    F[0] = 1;
    F[1] = 1;
} /* END InitRos2 */


/*----------------------------------------------------------------------------*/
/* A three-stage L-stable method of order 3                                   */
/*                                                                            */
/* E. Harier and G. Wanner, "Solving Ordenary Differential Equations II:      */
/* stiff and differential-algebraic problems." Computational Mathematics,     */
/* Springer-Verlag (1996)                                                     */
/*                                                                            */
/* @param[out] name     Method name                                           */
/* @param[out] nStage   Number of method stages                               */
/* @param[out] invLoEst One divided by the estimation of local order          */
/* @param[out] M        Coefficients for new step solution                    */
/* @param[out] E        Coefficients for error estimation                     */
/* @param[out] A        Lower triangular coefficient matrix                   */
/* @param[out] C        Lower triangular coefficient matrix                   */
/* @param[out] alpha    Y at stage i is approx. Y(T + H*Alpha_i)              */
/* @param[out] gamma    Stage i Gamma = sum(gamma[j])                         */
/* @param[out] F        Function evaluation flags                             */
/*----------------------------------------------------------------------------*/
void InitRos3(char ** name, int * nStage, float * invLoEst, float M[], float E[],
        float A[], float C[], float alpha[], float gamma[], char F[])
{
    /* Method name */
    *name = "Ros3";

    /* Number of stages */
    *nStage = 3;

    /* Inverse estimation of local order: 1/3 */
    *invLoEst = 0.3333333333333333;

    /* Coefficients for new step solution */
    M[0] = 1.0;
    M[1] = 6.1697947043828245592553615689730;
    M[2] = -0.4277225654321857332623837380651;

    /* Coefficients for error estimation */
    E[0] = 0.5;
    E[1] = -2.9079558716805469821718236208017;
    E[2] = 0.2235406989781156962736090927619;

    /* Lower triangular coefficient matrix A */
    A[0] = 1.0;
    A[1] = 1.0;
    A[2] = 0.0;

    /* Lower triangular coefficient matrix C */
    C[0] = -1.0156171083877702091975600115545;
    C[1] = 4.0759956452537699824805835358067;
    C[2] = 9.2076794298330791242156818474003;

    /* Two function evaluations */
    F[0] = 1;
    F[1] = 1;
    F[2] = 0;

    /* Y_stage_i ~ Y( T + H*Alpha_i ) */
    alpha[0] = 0.0;
    alpha[1] = 0.43586652150845899941601945119356;
    alpha[2] = 0.43586652150845899941601945119356;

    /* Gamma_i = \sum_j  gamma_{i,j}  */
    gamma[0] = 0.43586652150845899941601945119356;
    gamma[1] = 0.24291996454816804366592249683314;
    gamma[2] = 2.1851380027664058511513169485832;

} /* END InitRos3 */


/*----------------------------------------------------------------------------*/
/* A four-stage L-stable method of order 4                                    */
/*                                                                            */
/* E. Harier and G. Wanner, "Solving Ordenary Differential Equations II:      */
/* stiff and differential-algebraic problems." Computational Mathematics,     */
/* Springer-Verlag (1996)                                                     */
/*                                                                            */
/* @param[out] name     Method name                                           */
/* @param[out] nStage   Number of method stages                               */
/* @param[out] invLoEst One divided by the estimation of local order          */
/* @param[out] M        Coefficients for new step solution                    */
/* @param[out] E        Coefficients for error estimation                     */
/* @param[out] A        Lower triangular coefficient matrix                   */
/* @param[out] C        Lower triangular coefficient matrix                   */
/* @param[out] alpha    Y at stage i is approx. Y(T + H*Alpha_i)              */
/* @param[out] gamma    Stage i Gamma = sum(gamma[j])                         */
/* @param[out] F        Function evaluation flags                             */
/*----------------------------------------------------------------------------*/
void InitRos4(char ** name, int * nStage, float * invLoEst, float M[], float E[],
        float A[], float C[], float alpha[], float gamma[], char F[])
{
    /* Method name */
    *name = "Ros4";

    /* Number of stages */
    *nStage = 4;

    /* Inverse estimation of local order: 1/4 */
    *invLoEst = 0.25;

    /* Coefficients for new step solution */
    M[0] = 2.255570073418735;
    M[1] = 0.2870493262186792;
    M[2] = 0.4353179431840180;
    M[3] = 1.093502252409163;

    /* Coefficients for error estimation */
    E[0] = -0.2815431932141155;
    E[1] = -0.07276199124938920;
    E[2] = -0.1082196201495311;
    E[3] = -1.093502252409163;

    /* Lower triangular coefficient matrix A */
    A[0] = 2.0;
    A[1] = 1.867943637803922;
    A[2] = 0.2344449711399156;
    A[3] = 1.867943637803922;
    A[4] = 0.2344449711399156;
    A[5] = 0.0;

    /* Lower triangular coefficient matrix C */
    C[0] = -7.137615036412310;
    C[1] =  2.580708087951457;
    C[2] =  0.6515950076447975;
    C[3] = -2.137148994382534;
    C[4] = -0.3214669691237626;
    C[5] = -0.6949742501781779;

    /* Three function evaluations */
    F[0] = 1;
    F[1] = 1;
    F[2] = 1;
    F[3] = 0;

    /* Y_stage_i ~ Y( T + H*Alpha_i ) */
    alpha[0] = 0.0;
    alpha[1] = 1.145640000000000;
    alpha[2] = 0.6552168638155900;
    alpha[3] = 0.6552168638155900;

    /* Gamma_i = \sum_j  gamma_{i,j}  */
    gamma[0] = 0.5728200000000000;
    gamma[1] = -1.769193891319233;
    gamma[2] = 0.7592633437920482;
    gamma[3] = -0.1049021087100450;

} /* END InitRos4 */


/*----------------------------------------------------------------------------*/
/* A four-stage stiffly-stable method of order 4                              */
/*                                                                            */
/* E. Harier and G. Wanner, "Solving Ordenary Differential Equations II:      */
/* stiff and differential-algebraic problems." Computational Mathematics,     */
/* Springer-Verlag (1996)                                                     */
/*                                                                            */
/* @param[out] name     Method name                                           */
/* @param[out] nStage   Number of method stages                               */
/* @param[out] invLoEst One divided by the estimation of local order          */
/* @param[out] M        Coefficients for new step solution                    */
/* @param[out] E        Coefficients for error estimation                     */
/* @param[out] A        Lower triangular coefficient matrix                   */
/* @param[out] C        Lower triangular coefficient matrix                   */
/* @param[out] alpha    Y at stage i is approx. Y(T + H*Alpha_i)              */
/* @param[out] gamma    Stage i Gamma = sum(gamma[j])                         */
/* @param[out] F        Function evaluation flags                             */
/*----------------------------------------------------------------------------*/
void InitRodas3(char ** name, int * nStage, float * invLoEst, float M[], float E[],
        float A[], float C[], float alpha[], float gamma[], char F[])
{
    /* Method name */
    *name = "Rodas3";

    /* Number of stages */
    *nStage = 4;

    /* Inverse estimation of local order: 1/3 */
    *invLoEst = 0.3333333333333333;

    /* Coefficients for new step solution */
    M[0] = 2.0;
    M[1] = 0.0;
    M[2] = 1.0;
    M[3] = 1.0;

    /* Coefficients for error estimation */
    E[0] = 0.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 1.0;

    /* Lower triangular coefficient matrix A */
    A[0] = 0.0;
    A[1] = 2.0;
    A[2] = 0.0;
    A[3] = 2.0;
    A[4] = 0.0;
    A[5] = 1.0;

    /* Lower triangular coefficient matrix C */
    C[0] = 4.0;
    C[1] = 1.0;
    C[2] = -1.0;
    C[3] = 1.0;
    C[4] = -1.0;
    C[5] = -2.66666666666667;

    /* Three function evaluations */
    F[0] = 1;
    F[1] = 0;
    F[2] = 1;
    F[3] = 1;

    /* Y_stage_i ~ Y( T + H*Alpha_i ) */
    alpha[0] = 0.0;
    alpha[1] = 0.0;
    alpha[2] = 1.0;
    alpha[3] = 1.0;

    /* Gamma_i = \sum_j  gamma_{i,j}  */
    gamma[0] = 0.5;
    gamma[1] = 1.5;
    gamma[2] = 0.0;
    gamma[3] = 0.0;

} /* END InitRodas3 */


/*----------------------------------------------------------------------------*/
/* A six-stage stiffly-stable method of order 4                               */
/*                                                                            */
/* E. Harier and G. Wanner, "Solving Ordenary Differential Equations II:      */
/* stiff and differential-algebraic problems." Computational Mathematics,     */
/* Springer-Verlag (1996)                                                     */
/*                                                                            */
/* @param[out] name     Method name                                           */
/* @param[out] nStage   Number of method stages                               */
/* @param[out] invLoEst One divided by the estimation of local order          */
/* @param[out] M        Coefficients for new step solution                    */
/* @param[out] E        Coefficients for error estimation                     */
/* @param[out] A        Lower triangular coefficient matrix                   */
/* @param[out] C        Lower triangular coefficient matrix                   */
/* @param[out] alpha    Y at stage i is approx. Y(T + H*Alpha_i)              */
/* @param[out] gamma    Stage i Gamma = sum(gamma[j])                         */
/* @param[out] F        Function evaluation flags                             */
/*----------------------------------------------------------------------------*/
void InitRodas4(char ** name, int * nStage, float * invLoEst, float M[], float E[],
        float A[], float C[], float alpha[], float gamma[], char F[])
{
    /* Method name */
    *name = "Rodas4";

    /* Number of stages */
    *nStage = 6;

    /* Inverse estimation of local order: 1/4 */
    *invLoEst = 0.25;

    /* Coefficients for new step solution */
    M[0] = 1.544000000000000;
    M[1] = 6.019134481288629;
    M[2] = 12.53708332932087;
    M[3] = -0.6878860361058950;
    M[4] = 1.0;
    M[5] = 1.0;

    /* Coefficients for error estimation */
    E[0] = 0.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 0.0;
    E[4] = 0.0;
    E[5] = 1.0;

    /* Lower triangular coefficient matrix A */
    A[0] = 1.544000000000000;
    A[1] = 0.9466785280815826;
    A[2] = 0.2557011698983284;
    A[3] = 3.314825187068521;
    A[4] = 2.896124015972201;
    A[5] = 0.9986419139977817;
    A[6] = 1.221224509226641;
    A[7] = 6.019134481288629;
    A[8] = 12.53708332932087;
    A[9] = -0.6878860361058950;
    A[10] = 1.221224509226641;
    A[11] = 6.019134481288629;
    A[12] = 12.53708332932087;
    A[13] = -0.6878860361058950;
    A[14] = 1.0;

    /* Lower triangular coefficient matrix C */
    C[0] = -5.668800000000000;
    C[1] = -2.430093356833875;
    C[2] = -0.2063599157091915;
    C[3] = -0.1073529058151375;
    C[4] = -9.594562251023355;
    C[5] = -20.47028614809616;
    C[6] = 7.496443313967647;
    C[7] = -10.24680431464352;
    C[8] = -33.99990352819905;
    C[9] = 11.70890893206160;
    C[10] = 8.083246795921522;
    C[11] = -7.981132988064893;
    C[12] = -31.52159432874371;
    C[13] = 16.31930543123136;
    C[14] = -6.058818238834054;

    /* Six function evaluations */
    F[0] = 1;
    F[1] = 1;
    F[2] = 1;
    F[3] = 1;
    F[4] = 1;
    F[5] = 1;

    /* Y_stage_i ~ Y( T + H*Alpha_i ) */
    alpha[0] = 0.000;
    alpha[1] = 0.386;
    alpha[2] = 0.210;
    alpha[3] = 0.630;
    alpha[4] = 1.000;
    alpha[5] = 1.000;

    /* Gamma_i = \sum_j  gamma_{i,j}  */
    gamma[0] = 0.2500000000000000;
    gamma[1] = -0.1043000000000000;
    gamma[2] = 0.1035000000000000;
    gamma[3] = -0.03620000000000023;
    gamma[4] = 0.0;
    gamma[5] = 0.0;

} /* END InitRodas4 */


/*----------------------------------------------------------------------------*/
/* Calculates the left hand side matrix for Rosenbrock stage calculation      */
/*                                                                            */
/* @param[in]  diag     Value to add to diagonal elements                     */
/* @param[in]  jac      The Jacobian                                          */
/* @param[out] slhs     Left had side matrix for Rosenbrock stage calculation */
/*----------------------------------------------------------------------------*/
void RosenStageLHS(float diag, float * jac, float * slhs)
{
    int blk;
    int rem;
    float * slhs0 = slhs;
    
    blk = JAC_LU_NZ >> 3;
    rem = JAC_LU_NZ & 7;
    
    /* Copy and negate the Jacobian: slhs <== -jac */
    while(blk) {
        slhs[0] = -jac[0];
        slhs[1] = -jac[1];
        slhs[2] = -jac[2];
        slhs[3] = -jac[3];
        slhs[4] = -jac[4];
        slhs[5] = -jac[5];
        slhs[6] = -jac[6];
        slhs[7] = -jac[7];
        --blk;
        slhs += 8;
        jac += 8;
    }
    switch(rem) {
    case 7: slhs[6] = -jac[6];
    case 6: slhs[5] = -jac[5];
    case 5: slhs[4] = -jac[4];
    case 4: slhs[3] = -jac[3];
    case 3: slhs[2] = -jac[2];
    case 2: slhs[1] = -jac[1];
    case 1: slhs[0] = -jac[0];
    }
    
  /* Adjust diagonal elements */
  slhs0[0] = diag + slhs0[0];
  slhs0[3] = diag + slhs0[3];
  slhs0[8] = diag + slhs0[8];
  slhs0[13] = diag + slhs0[13];
  slhs0[18] = diag + slhs0[18];

} /* END RosenStageLHS */


/*----------------------------------------------------------------------------*/
/* Returns the "scaled norm" of the error vector                              */
/*----------------------------------------------------------------------------*/
float RosenErrNorm(float var[], float newY[], float errY[],
        size_t nTol, float const abstol[], float const reltol[])
{
    double err;
    float scale;
    float maxY;
    int i;

    err = ZERO;
    if (nTol > 1) {
        for (i=0; i<NVAR; ++i) {
            maxY = fmaxf(fabsf(var[i]),fabsf(newY[i]));
            scale = abstol[i] + reltol[i] * maxY;
            err = err + (errY[i] * errY[i]) / (scale * scale);
        }
    } else {
        for (i=0; i<NVAR; ++i) {
            maxY = fmaxf(fabsf(var[i]),fabsf(newY[i]));
            scale = abstol[0] + reltol[0] * maxY;
            err = err + (errY[i] * errY[i]) / (scale * scale);
        }
    }

    return (float)sqrt(err/(double)NVAR);

} /* END RosenErrNorm */



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
    int idata[20], float rdata[20])
{
    /* .................... Rosenbrock method parameters .................... */

    char * name;        /* Method name */
    int nStage;         /* Number of stages, from 2 to 6 */
    float invLoEst;     /* Inverse local order estimation */
    float M[6];         /* New step solution coefficients */
    float E[6];         /* Error estimation coefficients */
    float alpha[6];     /* Y_stage_i ~ Y( T + H*alpha_i ) */
    float gamma[6];     /* Gamma_i = \sum_j gamma_{i,j} */

    /* Coefficient matrices A and C are strictly lower triangular.
     * The subdiagonal elements are stored in row-wise order:
     * A(2,1)=A[0], A(3,1)=A[1], A(3,2)=A[2], etc. */
    float A[15];
    float C[15];

    /* F[i] == 0: stage i will re-use the function evaluation from stage i-1
     * F[i] != 0: stage i will evaluate function */
    char F[6];

    /* .................... Integration parameters .................... */

    float spanT;        /* Integration time span (positive value) */
    int autonomous;     /* idata[0]: Zero if F = F(t,y) */
    int nTol;           /* idata[1]: Length of the tolerance vectors, 1 = scalar*/
    int stepMax;        /* idata[2]: Maximum permitted steps */
    float minH;         /* rdata[0]: Integration step size lower bound */
    float maxH;         /* rdata[1]: Integration step size upper bound */
    float startH;       /* rdata[2]: Starting integration step size */
    float minFact;      /* rdata[3]: Lower bound on step decrease factor */
    float maxFact;      /* rdata[4]: Upper bound on step increase factor */
    float rejectFact;   /* rdata[5]: Step decrease factor after step rejection */
    float safeFact;     /* rdata[6]: Safety factor in computation of new step size */

    /* .................... Local variables .................... */

    float * K = 0;          /* Stage solution vectors */
    float * newVar = 0;     /* Variable concentrations after successful solve */
    float * errVar = 0;     /* Error in newVar */
    float * fcn0 = 0;       /* Function at time tstart */
    float * fcn = 0;        /* Function at time T */
    float * dFdT = 0;       /* Partial derivative of the function w.r.t T */
    float * rct = 0;        /* Reaction rates at time T */
    float * jac0 = 0;       /* Jacobian at time tstart */
    float * slhs = 0;       /* Stage computation left hand side matrix */

    int dir;            /* +1 if time advances positively, -1 otherwise */
    float T;            /* Model time */
    float H;            /* Timestep */
    float newH;         /* Updated timestep */
    float errNorm;      /* The scaled norm of the error vector */

    int rejectH;        /* Number of consecutive time step rejections */
    int i, j;           /* Iterators */

    int nFun = 0;       /* Number of function evaluations */
    int nJac = 0;       /* Number of Jacobian evaluations */
    int nStp = 0;       /* Number of solver steps */
    int nAcc = 0;       /* Number of accepted steps */
    int nRej = 0;       /* Number of rejected steps */
    int nDec = 0;       /* Number of matrix decompositions */
    int nSol = 0;       /* Number of Ax=b solves */
    int nSng = 0;       /* Number of singular decomposition results */

    /* Macro to clean up and abort the integrator */
    #define ABORT(code, fmt, ...) { \
        printf("Kppa: %s: T=%g, H=%g: " fmt, name, T, H, ##__VA_ARGS__); \
        idata[19] = code; \
        goto end; \
    }

    /* ................ Initialize the Rosenbrock method ................ */

    name = "Unknown";
    switch (idata[3]) {
    case 0:
        InitRos4(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    case 1:
        InitRos2(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    case 2:
        InitRos3(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    case 3:
        InitRos4(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    case 4:
        InitRodas3(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    case 5:
        InitRodas4(&name, &nStage, &invLoEst, M, E, A, C, alpha, gamma, F);
        break;
    default:
        fprintf(stderr, "Kppa: Unknown method: %d\n", idata[3]);
        idata[19] = -3;
        return;
    }

    /* ................... Initialize local variables ................... */

    /* Initialize step rejection counter */
    rejectH = 0;

    /* Initialize time */
    dir = (tend >= tstart ? +1 : -1);
    spanT = dir * (tend - tstart);
    T = tstart;
    H = spanT;

    /* Determine if F depends on time */
    autonomous = (idata[0] != 0);

    /* Scalar tolerances limits the tolerance vectors to the first element. */
    nTol = idata[1] ? 1 : NVAR;
    
    /* Initialize error norm */
    errNorm = (float)1.0 / FLT_EPSILON;

    /* Maximum number of steps before the method aborts */
    stepMax = idata[2] ? idata[2] : 100000;
    if (stepMax < 0)
        ABORT(-3, "Invalid maximum steps: %d\n", stepMax);

    /* Check tolerance vectors */
    if(idata[4]) {
        for (i=0; i<nTol; i++) {
            if (abstol[i] <= ZERO)
                ABORT(-3, "Unreasonable tolerance: abstol[%d]=%g\n", i, abstol[i]);
            if (reltol[i] <= (10.0 * FLT_EPSILON) || reltol[i] >= ONE)
                ABORT(-3, "Unreasonable tolerance: reltol[%d]=%g\n", i, reltol[i]);
        }
    }

    /* Lower bound on the step size: (positive value) */
    minH = rdata[0];
    if (minH < ZERO)
        ABORT(-3, "Invalid step size lower bound: %g\n", minH);

    /* Upper bound on the step size: (positive value) */
    maxH = rdata[1] ? fminf(fabsf(rdata[1]), spanT) : spanT;
    if (maxH < ZERO)
        ABORT(-3, "Invalid step size upper bound: %g\n", maxH);

    /*  Starting step size: (positive value) */
    startH = rdata[2] ? fminf(fabsf(rdata[2]), spanT) : fmaxf(minH,MIN_DELT);
    if (startH < ZERO)
        ABORT(-3, "Invalid starting step size: %g\n", startH);

    /* Lower bound on step decrease factor */
    minFact = rdata[3] ? rdata[3] : 0.2;
    if (minFact < ZERO)
        ABORT(-3, "Invalid lower bound on step decrease factor: %g\n", minFact);

    /* Upper bound on step increase factor */
    maxFact = rdata[4] ? rdata[4] : 6.0;
    if (maxFact < minFact)
        ABORT(-3, "Invalid upper bound on step increase factor: %g\n", maxFact);

    /* Step decrease factor after step rejection */
    rejectFact = rdata[5] ? rdata[5] : 0.1;
    if (rejectFact < ZERO)
        ABORT(-3, "Invalid step decrease factor for rejected step: %g\n", rejectFact);

    /* Safety factor in the computation of new step size */
    safeFact = rdata[6] ? rdata[6] : 0.9;
    if (safeFact < ZERO)
        ABORT(-3, "Invalid new step safety factor: %g\n", safeFact);

    /* Adjust timestep according to user-specified limits */
    H = fminf(startH, maxH);
    if (fabsf(H) < 10 * FLT_EPSILON)
        H = MIN_DELT;
        
    /* Allocate memory */
    K = (float*)malloc(nStage*NVAR*sizeof(float));
    newVar = (float*)malloc(NVAR*sizeof(float));
    errVar = (float*)malloc(NVAR*sizeof(float));
    fcn0 = (float*)malloc(NVAR*sizeof(float));
    fcn  = (float*)malloc(NVAR*sizeof(float));
    dFdT = (float*)malloc(NVAR*sizeof(float));
    rct  = (float*)malloc(NREACT*sizeof(float));
    jac0 = (float*)malloc(JAC_LU_NZ*sizeof(float));
    slhs = (float*)malloc(JAC_LU_NZ*sizeof(float));

    /* ............................ Integrate ............................ */

    while(fabsf(tend - T) > FLT_EPSILON) {   /* Time integration loop */

        /* Check step count */
        if (nStp > stepMax)
            ABORT(-6, "Too many integration steps: stepMax=%d\n", stepMax);

        /* Check timestep size */
        if ((T + 0.1*H == T) || (H <= FLT_EPSILON))
            ABORT(-7, "Step size too small (T + H/10 = T) or H < eps\n");

        /* Update timestep */
        H = fminf(H,fabsf(tend-T));

        /* Compute reaction rates at the current time */
        Rates(T, idx, rct);

        /* Compute the function at the current time */
        Fun(var, fix, rct, fcn0);
        ++nFun;

        /* Compute the Jacobian at the current time */
        Jac(var, fix, rct, jac0);
        ++nJac;

        /* Compute the function derivative with respect to time */
        if (!autonomous) {
            float delta = sqrt(FLT_EPSILON) * fmaxf(MIN_DELT, fabsf(T));
            Rates(T+delta, idx, rct);
            Fun(var, fix, rct, fcn);
            ++nFun;
            WYMXDA(NVAR, fcn0, fcn, delta, dFdT);
        }

        /* Repeat step calculation until step accepted  */
        do {
            int singRow = 0;
            int decomps = 1;

            /* Prepare the LHS matrix for stage calculations */
            RosenStageLHS(1.0/(dir*H*gamma[0]), jac0, slhs);

            /* LU decompose stage LHS matrix */
            singRow = Decomp(slhs);
            ++nDec;

            /* If the decomposition failed, half the timestep and try again */
            while(singRow) {
                printf("Kppa: %s: LU decomposition singular on row %d\n", name, singRow-1);
                ++nSng;

                /* Reduce step size */
                H *= HALF;

                /* Abort after eight failed decompositions */
                if (decomps > 8 || H == ZERO)
                    ABORT(-8, "Matrix is repeatedly singular\n");

                /* Build new stage LHS with reduced time step */
                RosenStageLHS(1.0/(dir*H*gamma[0]), jac0, slhs);

                /* LU decompose stage LHS matrix */
                singRow = Decomp(slhs);
                ++nDec;
                ++decomps;
            }

            /* Compute stage 0 using the previously-computed function */
            WCOPY(NVAR, fcn0, fcn);
            WCOPY(NVAR, fcn0, K);
            if (!autonomous && gamma[0]) {
                WAXPY(NVAR, dir*H*gamma[0], dFdT, K);
            }

            /* Solve stage 0 */
            Solve(slhs, K);
            nSol++;

            /* Compute the remaining stages  */
            for (i=1; i<nStage; ++i) {
                float * Kstage = K + i*NVAR;

                if (F[i]) {
                    float tau = T + alpha[i] * dir * H;
                    /* Apply coefficient matrix A */
                    WCOPY(NVAR, var, newVar);
                    for (j=0; j<i; ++j) {
                        WAXPY(NVAR, A[i*(i-1)/2+j], K+j*NVAR, newVar);
                    }
                    /* Update reaction rates, if necessary */
                    if(!autonomous) {
                        Rates(tau, idx, rct);
                    }
                    /* Evaluate the function */
                    Fun(newVar, fix, rct, fcn);
                    ++nFun;
                }

                /* Apply coefficient matrix C */
                WCOPY(NVAR, fcn, Kstage);
                for (j=0; j<i; ++j) {
                    WAXPY(NVAR, C[i*(i-1)/2+j]/(dir*H), K+j*NVAR, Kstage);
                }

                if (!autonomous && gamma[i]) {
                    WAXPY(NVAR, dir*H*gamma[i], dFdT, Kstage);
                }

                /* Solve stage i */
                Solve(slhs, Kstage);
                nSol++;
            }

            /* Compute the new solution */
            WCOPY(NVAR, var, newVar);
            for (j=0; j<nStage; j++) {
                WAXPY(NVAR, M[j], K+j*NVAR, newVar);
            }

            /* Estimate error */
            WSCAL(NVAR, ZERO, errVar);
            for (j=0; j<nStage; j++) {
                WAXPY(NVAR, E[j], K+j*NVAR, errVar);
            }
            
            /* Calculate scaled norm of the error vector */
            errNorm = RosenErrNorm(var, newVar, errVar, nTol, abstol, reltol);
            rdata[12] = errNorm;
            if(isinf(errNorm)) {
                ABORT(-10, "Error norm in cell %ld is Inf\n", idx);
            } else if(isnan(errNorm)) {
                ABORT(-10, "Error norm in cell %ld is NaN\n", idx);
            }

            /* Calculate a new step size: minFact <= newH/H <= maxFact */
            newH = H * fminf(maxFact,fmaxf(minFact,safeFact/pow(errNorm,invLoEst)));
            ++nStp;

            /* Decide to accept or reject step  */
            if (errNorm <= ONE || H <= minH) {
                /* Step accepted */
                ++nAcc;
                WCOPY(NVAR, newVar, var);
                T += dir * H;
                /* Adjust step size */
                newH = fmaxf(minH,fminf(newH,maxH));
                if(rejectH) {
                    newH = fminf(newH,H);
                }
                rejectH = 0;
                H = newH;
                /* Return to time loop */
                break;
            } else {
                /* Step rejected */
                ++nRej;
                if(rejectH > 1) {
                    newH = H * rejectFact;
                }
                ++rejectH;
                H = newH;
                /* Continue step calculation */
                continue;
            }

        } while(1); /* Step calculation */

    } /* Time loop */

    /* ...................... Exit integrator ...................... */
    
    /* Set exit status */
    idata[19] = 0;

end:

    /* Deallocate memory */
    free((void*)K);
    free((void*)newVar);
    free((void*)errVar);
    free((void*)fcn0);
    free((void*)fcn);
    free((void*)rct);
    free((void*)dFdT);
    free((void*)jac0);
    free((void*)slhs);

    /* Collect statistics */
    idata[10] = nFun;
    idata[11] = nJac;
    idata[12] = nStp;
    idata[13] = nAcc;
    idata[14] = nRej;
    idata[15] = nDec;
    idata[16] = nSol;
    idata[17] = nSng;
    /* Record exit time and last step size */
    rdata[10] = T;
    rdata[11] = H;
    rdata[12] = errNorm;

}/* END Integrate */



/*------------------------- END small_rosenbrock.h END ------------------------*/
