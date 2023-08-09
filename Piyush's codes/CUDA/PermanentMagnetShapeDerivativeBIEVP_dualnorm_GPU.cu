#include <cuda_runtime.h>
#include <stdio.h>
#include <cuda.h>
#include <math.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

__device__ double atomicAdd(double *address, double val)
{
    unsigned long long int *address_as_ull =
        (unsigned long long int *)address;
    unsigned long long int old = *address_as_ull, assumed;

    do
    {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

__device__ double SLKernel(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    return 1. / (4 * M_PI) / YmX.norm();
}

__device__ Eigen::Vector3d Vel(int a, int b, int c, int alpha, Eigen::Vector3d X)
{
    Eigen::Vector3d out = Eigen::Vector3d::Zero(3);
    out(alpha) = cos(a * X(0)) * cos(b * X(1)) * cos(c * X(2));
    return out;
}

__device__ Eigen::Matrix3d DVel(int a, int b, int c, int alpha, Eigen::Vector3d X)
{
    Eigen::Matrix3d out = Eigen::Matrix3d::Zero(3, 3);
    out.row(alpha) = -Eigen::Vector3d(a * sin(a * X(0)) * cos(b * X(1)) * cos(c * X(2)),
                                      b * cos(a * X(0)) * sin(b * X(1)) * cos(c * X(2)),
                                      c * cos(a * X(0)) * cos(b * X(1)) * sin(c * X(2)));
    return out;
}

__device__ double KernelA1(int a, int b, int c, int alpha, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    return YmX.dot(Vel(a, b, c, alpha, X) - Vel(a, b, c, alpha, Y)) / (4 * M_PI) / (znorm * znorm * znorm);
}

__device__ double KernelA2(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    return 1. / (4 * M_PI) / YmX.norm();
}

__device__ Eigen::Vector3d KernelC1(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    double znorm = YmX.norm();
    return YmX / (4 * M_PI) / (znorm * znorm * znorm);
}

__device__ Eigen::Vector3d KernelC3(int a, int b, int c, int alpha, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    double znorm = YmX.norm();
    return -3. / (4 * M_PI) * YmX * YmX.dot(Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) / (znorm * znorm * znorm * znorm * znorm) //
           + 1. / (4 * M_PI) * (Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) / (znorm * znorm * znorm);
}

__device__ Eigen::Vector3d DLKernel(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    double znorm = YmX.norm();
    return -YmX / (4 * M_PI) / (znorm * znorm * znorm);
}

// Returns Intersection, DiffI, DiffJ
__device__ void IntersectionDiff(int *EltI, int *EltJ, int intersection[], int diffI[], int diffJ[])
{
    bool EltITracker[] = {false, false, false};
    bool EltJTracker[] = {false, false, false};

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (EltI[i] == EltJ[j])
            {
                EltITracker[i] = true;
                EltJTracker[j] = true;
            }
        }
    }

    int interidx = 0, diffiidx = 0, diffjidx = 0;

    for (int i = 0; i < 3; ++i)
    {
        if (EltITracker[i])
            intersection[interidx++] = EltI[i];
        else
            diffI[diffiidx++] = EltI[i];

        if (!EltJTracker[i])
            diffJ[diffjidx++] = EltJ[i];
    }

    /* for (int i = 0; i < 3; ++i)
    {
        if (EltJTracker[i])
            intersection[interidx++] = EltJ[i];
        else
            diffJ[diffjidx++] = EltJ[i];

        if (!EltITracker[i])
            diffI[diffiidx++] = EltI[i];
    } */
}

__global__ void computeShapeDerivative(int TrialDim, int TestDim, int NTriangles, int NVertices, int NInteractions,
                                       int NThreads, const unsigned short int *I, const unsigned short int *J, const int *relation,
                                       const double *W0, const double *X0, int Nq0,
                                       const double *W1, const double *X1, int Nq1,
                                       const double *W2, const double *X2, int Nq2,
                                       const double *W3, const double *X3, int Nq3,
                                       double *A1, double *A2, double *C1, double *C2, double *C3, double *N,
                                       double *red_remaining, double *red_l_M, double *blue_remaining, double *blue_l_M,
                                       double *MC1, double *MC2, double *MC3,
                                       double *shapeDerivative,
                                       const double *TdA, const double *TnA, const double *Mxn_coeffs,
                                       const double mu, const double mu0,
                                       const int *Elements, const double *Vertices, const double *Normals, const double *Areas,
                                       const int *Elt2DofTest, const int *Elt2DofTrial,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest,
                                       const int *abc_alpha, int Nabc_alpha,
                                       const int *permII, const int *permJJ)
/* double *testout, double *testABCi, double *testABCj,
int *orig_elti, int *orig_eltj, int *modif_elti, int *modif_eltj) */
{
    // shared memory for writing shape derivatives
    // Need to specify the size beforehand, hard coding Nabc_alpha = 81, blockDim.x = 32
    __shared__ double localShapeDerivatives[81 * 32];
    // Initialization of shared memory
    for (int l = 0; l < Nabc_alpha; ++l)
    {
        localShapeDerivatives[l + threadIdx.x * Nabc_alpha] = 0;
    }
    int ThreadID = blockIdx.x * blockDim.x + threadIdx.x;

    //*shapeDerivative = 3.145;

    // Size of the matrix
    // int NInteractions = TrialDim * TestDim;

    // Number of element interaction
    // int NInteractions = NTriangles * NTriangles;

    // Matrix size is NTriangles^2. Each entry is assigned to a thread
    // InteractionsPerThread gives the max no.of element interactions assigned
    // to a thread
    int InteractionsPerThread = ceil(double(NInteractions) / double(NThreads));

    /* if (blockIdx.x == 0 && threadIdx.x == 0)
    {
        printf("Number of blocks: %d \n ", gridDim.x);
        printf("Threads per block: %d \n ", blockDim.x);
        printf("Total interactions: %d , Interactions per thread: %d \n ", NInteractions, InteractionsPerThread);
    } */

    // Looping over all assigned interactions
    for (int idx = 0; idx < InteractionsPerThread; ++idx)
    {
        /* if (blockIdx.x == 0 && threadIdx.x == 0)
        {
            printf("In block 0 thread 0 computing interaction no. : %d \n ", idx);
        } */
        // The interaction number
        // int InteractionIdx = ThreadID + NThreads * idx;
        int InteractionIdx = ThreadID * InteractionsPerThread + idx;

        // Preparing variables
        Eigen::Vector3d Ai, Bi, Ci, Aj, Bj, Cj;
        Eigen::MatrixXd Ei(3, 2), Ej(3, 2);

        // Make sure that the last thread stays in limit
        if (InteractionIdx >= NInteractions)
            break;

        // The pair of panels
        int i = __ldg(&I[InteractionIdx]), j = __ldg(&J[InteractionIdx]);

        /* printf("Interaction  (%d, %d) \n", i, j); */

        double g_tau = 2 * Areas[i], g_t = 2 * Areas[j];

        Eigen::Vector3d normalx(__ldg(&Normals[3 * i]), __ldg(&Normals[3 * i + 1]), __ldg(&Normals[3 * i + 2]));
        Eigen::Vector3d normaly(__ldg(&Normals[3 * j]), __ldg(&Normals[3 * j + 1]), __ldg(&Normals[3 * j + 2]));

        int EltI[] = {__ldg(&Elements[3 * i]),
                      __ldg(&Elements[3 * i + 1]),
                      __ldg(&Elements[3 * i + 2])};

        int EltJ[] = {__ldg(&Elements[3 * j]),
                      __ldg(&Elements[3 * j + 1]),
                      __ldg(&Elements[3 * j + 2])};

        int origEltI[] = {EltI[0], EltI[1], EltI[2]};
        int origEltJ[] = {EltJ[0], EltJ[1], EltJ[2]};

        /* int DofsI[] = {Elt2DofTest[i], Elt2DofTest[i + NTriangles], Elt2DofTest[i + 2 * NTriangles]};
        int DofsJ[] = {Elt2DofTrial[j], Elt2DofTrial[j + NTriangles], Elt2DofTrial[j + 2 * NTriangles]}; */

        int DofsI[] = {__ldg(&Elt2DofTest[3 * i]), __ldg(&Elt2DofTest[3 * i + 1]), __ldg(&Elt2DofTest[3 * i + 2])};
        int DofsJ[] = {__ldg(&Elt2DofTrial[3 * j]), __ldg(&Elt2DofTrial[3 * j + 1]), __ldg(&Elt2DofTrial[3 * j + 2])};

        // Original permutation of elements
        int permI[] = {0, 1, 2};
        int permJ[] = {0, 1, 2};

        const double *Weights[4] = {W0, W1, W2, W3};
        const double *Points[4] = {X0, X1, X2, X3};
        const int NumPoints[4] = {Nq0, Nq1, Nq2, Nq3};

        const double *W = Weights[relation[InteractionIdx]];
        const double *X = Points[relation[InteractionIdx]];
        int NQudPts = NumPoints[relation[InteractionIdx]];

        /* if (relation[InteractionIdx] == 0) // No interaction
        {
            // Computing Quadrature
            W = W0;
            X = X0;
            NQudPts = Nq0;
        }
        else if (relation[InteractionIdx] == 1) // Common vertex
        {
            // Computing Quadrature
            W = W1;
            X = X1;
            NQudPts = Nq1;
        }
        else if (relation[InteractionIdx] == 2) // Common edge
        {
            // Computing Quadrature
            W = W2;
            X = X2;
            NQudPts = Nq2;
        }
        else // Identical panels, case 3
        {
            // Computing Quadrature
            W = W3;
            X = X3;
            NQudPts = Nq3;
        } */

        // EltI and EltJ changed according to permII and permJJ
        for (int k = 0; k < 3; ++k)
        {
            permI[k] = permII[3 * InteractionIdx + k];
            permJ[k] = permJJ[3 * InteractionIdx + k];
            EltI[k] = origEltI[permI[k]];
            EltJ[k] = origEltJ[permJ[k]];
        }

        // Vertices of element i
        Ai = Eigen::Vector3d(__ldg(&Vertices[3 * EltI[0]]), __ldg(&Vertices[3 * EltI[0] + 1]), __ldg(&Vertices[3 * EltI[0] + 2]));
        Bi = Eigen::Vector3d(__ldg(&Vertices[3 * EltI[1]]), __ldg(&Vertices[3 * EltI[1] + 1]), __ldg(&Vertices[3 * EltI[1] + 2]));
        Ci = Eigen::Vector3d(__ldg(&Vertices[3 * EltI[2]]), __ldg(&Vertices[3 * EltI[2] + 1]), __ldg(&Vertices[3 * EltI[2] + 2]));

        // Vertices of element j
        Aj = Eigen::Vector3d(__ldg(&Vertices[3 * EltJ[0]]), __ldg(&Vertices[3 * EltJ[0] + 1]), __ldg(&Vertices[3 * EltJ[0] + 2]));
        Bj = Eigen::Vector3d(__ldg(&Vertices[3 * EltJ[1]]), __ldg(&Vertices[3 * EltJ[1] + 1]), __ldg(&Vertices[3 * EltJ[1] + 2]));
        Cj = Eigen::Vector3d(__ldg(&Vertices[3 * EltJ[2]]), __ldg(&Vertices[3 * EltJ[2] + 1]), __ldg(&Vertices[3 * EltJ[2] + 2]));

        // Jacobian Matrices

        Ei.col(0) = Bi - Ai;
        Ei.col(1) = Ci - Ai;

        Ej.col(0) = Bj - Aj;
        Ej.col(1) = Cj - Aj;

        Eigen::MatrixXd EtEi = Ei.transpose() * Ei;
        Eigen::MatrixXd EtEj = Ej.transpose() * Ej;
        double deti = EtEi(0, 0) * EtEi(1, 1) - EtEi(0, 1) * EtEi(1, 0);
        double detj = EtEj(0, 0) * EtEj(1, 1) - EtEj(0, 1) * EtEj(1, 0);

        Eigen::MatrixXd Dxyi = -EtEi;
        Eigen::MatrixXd Dxyj = -EtEj;

        Dxyi(0, 0) = EtEi(1, 1);
        Dxyi(1, 1) = EtEi(0, 0);

        Dxyj(0, 0) = EtEj(1, 1);
        Dxyj(1, 1) = EtEj(0, 0);

        Dxyi /= deti;
        Dxyj /= detj;

        Eigen::MatrixXd DCVi = Ei * Dxyi, DCVj = Ej * Dxyj;

        // RWG X RWG
        /* for (int ii = 0; ii < 3; ++ii)
        {
            int iip1 = (permI[ii] + 1) % 3;
            int iip2 = (iip1 + 1) % 3;

            double fluxI = origEltI[iip1] < origEltI[iip2] ? 1. : -1.;

            for (int jj = 0; jj < 3; ++jj)
            {
                int jjp1 = (permJ[jj] + 1) % 3;
                int jjp2 = (jjp1 + 1) % 3;

                double fluxJ = origEltJ[jjp1] < origEltJ[jjp2] ? 1. : -1.;

                for (int fieldIdx = 0; fieldIdx < Nabc_alpha; ++fieldIdx)
                {
                    Eigen::Matrix3d LocalMatrix = Eigen::MatrixX3d::Zero(3, 3);
                    Eigen::Matrix3d LocalMatrixA1 = Eigen::MatrixX3d::Zero(3, 3);
                    Eigen::Matrix3d LocalMatrixA2 = Eigen::MatrixX3d::Zero(3, 3);
                    Eigen::Matrix3d LocalMatrixC1 = Eigen::MatrixX3d::Zero(3, 3);
                    Eigen::Matrix3d LocalMatrixC3 = Eigen::MatrixX3d::Zero(3, 3);
                    Eigen::Matrix3d LocalMatrixN = Eigen::MatrixX3d::Zero(3, 3);

                    // Computing the shape derivaive for the field specified by fieldIdx
                    for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                    {

                        double RWGX_ref_0 = __ldg(&X[4 * QudPt]) - ii % 2;
                        double RWGX_ref_1 = __ldg(&X[4 * QudPt + 1]) - ii / 2;

                        double RWGY_ref_0 = __ldg(&X[4 * QudPt + 2]) - jj % 2;
                        double RWGY_ref_1 = __ldg(&X[4 * QudPt + 3]) - jj / 2;

                        // RWG elements
                        Eigen::Vector3d Psix = fluxI * (Ei.col(0) * RWGX_ref_0 + Ei.col(1) * RWGX_ref_1);
                        Eigen::Vector3d Psiy = fluxJ * (Ej.col(0) * RWGY_ref_0 + Ej.col(1) * RWGY_ref_1);

                        Eigen::Vector3d chi_tau = Ai + Ei.col(0) * __ldg(&X[4 * QudPt]) + Ei.col(1) * __ldg(&X[4 * QudPt + 1]);
                        Eigen::Vector3d chi_t = Aj + Ej.col(0) * __ldg(&X[4 * QudPt + 2]) + Ej.col(1) * __ldg(&X[4 * QudPt + 3]);

                        // A1 with RWG . RWG
                        LocalMatrixA1(ii, jj) += __ldg(&W[QudPt]) * KernelA1(__ldg(&abc_alpha[4 * fieldIdx]), __ldg(&abc_alpha[4 * fieldIdx + 1]), __ldg(&abc_alpha[4 * fieldIdx + 2]), __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau) * Psiy.dot(Psix);

                        // A2 with DVelRWG(y) . RWG(x)
                        LocalMatrixA2(ii, jj) += __ldg(&W[QudPt]) * KernelA2(chi_tau, chi_t, chi_t - chi_tau) * (DVel(__ldg(&abc_alpha[4 * fieldIdx]), __ldg(&abc_alpha[4 * fieldIdx + 1]), __ldg(&abc_alpha[4 * fieldIdx + 2]), __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_t) * Psiy).dot(Psix);

                        // C1 with DVelRWG(y) X RWG(X)
                        LocalMatrixC1(ii, jj) += __ldg(&W[QudPt]) * (KernelC1(chi_tau, chi_t, chi_t - chi_tau).cross(DVel(__ldg(&abc_alpha[4 * fieldIdx]), __ldg(&abc_alpha[4 * fieldIdx + 1]), __ldg(&abc_alpha[4 * fieldIdx + 2]), __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_t) * Psiy)).dot(Psix);

                        // C3 with RWG X RWG
                        LocalMatrixC3(ii, jj) += __ldg(&W[QudPt]) * (KernelC3(__ldg(&abc_alpha[4 * fieldIdx]), __ldg(&abc_alpha[4 * fieldIdx + 1]), __ldg(&abc_alpha[4 * fieldIdx + 2]), __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau).cross(Psiy)).dot(Psix);

                        // N with RWG.div X RWG.div, kernelN = kernelA1
                        LocalMatrixN(ii, jj) += __ldg(&W[QudPt]) * KernelA1(__ldg(&abc_alpha[4 * fieldIdx]), __ldg(&abc_alpha[4 * fieldIdx + 1]), __ldg(&abc_alpha[4 * fieldIdx + 2]), __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau) * 4 * fluxI * fluxJ;
                    }

                    // Adding the shape derivative to the local shared storage
                    localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += 1 / (2 * mu0) * ((1 + mu / mu0) * (TnA[DofsI[permI[ii]]] * LocalMatrixA1(ii, jj) * TnA[DofsJ[permJ[jj]]]        // A1
                                                                                                                     + 2 * TnA[DofsI[permI[ii]]] * LocalMatrixA2(ii, jj) * TnA[DofsJ[permJ[jj]]]) // A2
                                                                                                   - 4 * (TnA[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * TdA[DofsJ[permJ[jj]]]                   // C1
                                                                                                          + TdA[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * TnA[DofsJ[permJ[jj]]]                 // C2
                                                                                                          + TnA[DofsI[permI[ii]]] * LocalMatrixC3(ii, jj) * TdA[DofsJ[permJ[jj]]])                // C3
                                                                                                   + (1 + mu0 / mu) * -TdA[DofsI[permI[ii]]] * LocalMatrixN(ii, jj) * TdA[DofsJ[permJ[jj]]])      // N
                                                                                  + mu / 2 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA1(ii, jj) * Mxn_coeffs[DofsJ[permJ[jj]]]                  // red_remaining
                                                                                  + mu * TnA[DofsI[permI[ii]]] * LocalMatrixA1(ii, jj) * Mxn_coeffs[DofsJ[permJ[jj]]]                             // red_l_M
                                                                                  + mu * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2(ii, jj) * Mxn_coeffs[DofsJ[permJ[jj]]]                      // blue_remaining
                                                                                  + mu * (Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2(ii, jj) * TnA[DofsJ[permJ[jj]]] +                          // blue_l_M
                                                                                          TnA[DofsI[permI[ii]]] * LocalMatrixA2(ii, jj) * Mxn_coeffs[DofsJ[permJ[jj]]]) +
                                                                                  -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * TdA[DofsJ[permJ[jj]]]    // MC1
                                                                                  + -mu0 * TdA[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * Mxn_coeffs[DofsJ[permJ[jj]]]  // MC2
                                                                                  + -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC3(ii, jj) * TdA[DofsJ[permJ[jj]]]; // MC3
                }
            }
        } */

        // Changed loop order, MUCH SLOWER!!!
        for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
        {
            Eigen::Vector3d chi_tau = Ai + Ei.col(0) * __ldg(&X[4 * QudPt]) + Ei.col(1) * __ldg(&X[4 * QudPt + 1]);
            Eigen::Vector3d chi_t = Aj + Ej.col(0) * __ldg(&X[4 * QudPt + 2]) + Ej.col(1) * __ldg(&X[4 * QudPt + 3]);

            // Looping over different fields
            for (int fieldIdx = 0; fieldIdx < Nabc_alpha; ++fieldIdx)
            {
                double localval = 0;
                for (int ii = 0; ii < 3; ++ii)
                {
                    int iip1 = (permI[ii] + 1) % 3;
                    int iip2 = (iip1 + 1) % 3;

                    double fluxI = origEltI[iip1] < origEltI[iip2] ? 1. : -1.;
                    for (int jj = 0; jj < 3; ++jj)
                    {
                        int jjp1 = (permJ[jj] + 1) % 3;
                        int jjp2 = (jjp1 + 1) % 3;

                        double fluxJ = origEltJ[jjp1] < origEltJ[jjp2] ? 1. : -1.;

                        double RWGX_ref_0 = __ldg(&X[4 * QudPt]) - ii % 2;
                        double RWGX_ref_1 = __ldg(&X[4 * QudPt + 1]) - ii / 2;

                        double RWGY_ref_0 = __ldg(&X[4 * QudPt + 2]) - jj % 2;
                        double RWGY_ref_1 = __ldg(&X[4 * QudPt + 3]) - jj / 2;

                        // RWG elements
                        Eigen::Vector3d Psix = fluxI * (Ei.col(0) * RWGX_ref_0 + Ei.col(1) * RWGX_ref_1);
                        Eigen::Vector3d Psiy = fluxJ * (Ej.col(0) * RWGY_ref_0 + Ej.col(1) * RWGY_ref_1);

                        // A1 with RWG . RWG
                        double LocalMatrixA1_ii_jj = __ldg(&W[QudPt]) //
                                                     * KernelA1(abc_alpha[4 * fieldIdx],
                                                                abc_alpha[4 * fieldIdx + 1],
                                                                abc_alpha[4 * fieldIdx + 2],
                                                                abc_alpha[4 * fieldIdx + 3],
                                                                chi_tau, chi_t, chi_t - chi_tau) *
                                                     Psiy.dot(Psix);

                        // A2 with DVelRWG(y) . RWG(x)
                        double LocalMatrixA2_ii_jj = __ldg(&W[QudPt])                                                                 //
                                                     * KernelA2(chi_tau, chi_t, chi_t - chi_tau) * (DVel(abc_alpha[4 * fieldIdx],     //
                                                                                                         abc_alpha[4 * fieldIdx + 1], //
                                                                                                         abc_alpha[4 * fieldIdx + 2], //
                                                                                                         abc_alpha[4 * fieldIdx + 3], //
                                                                                                         chi_t) *
                                                                                                    Psiy)
                                                                                                       .dot(Psix);

                        // C1 with DVelRWG(y) X RWG(X)
                        double LocalMatrixC1_ii_jj = __ldg(&W[QudPt]) *                                                                 //
                                                     (KernelC1(chi_tau, chi_t, chi_t - chi_tau).cross(DVel(abc_alpha[4 * fieldIdx],     //
                                                                                                           abc_alpha[4 * fieldIdx + 1], //
                                                                                                           abc_alpha[4 * fieldIdx + 2], //
                                                                                                           abc_alpha[4 * fieldIdx + 3], chi_t) *
                                                                                                      Psiy))
                                                         .dot(Psix);

                        // C3 with RWG X RWG
                        double LocalMatrixC3_ii_jj = __ldg(&W[QudPt])                         //
                                                     * (KernelC3(abc_alpha[4 * fieldIdx],     //
                                                                 abc_alpha[4 * fieldIdx + 1], //
                                                                 abc_alpha[4 * fieldIdx + 2], //
                                                                 abc_alpha[4 * fieldIdx + 3], //
                                                                 chi_tau, chi_t, chi_t - chi_tau)
                                                            .cross(Psiy))
                                                           .dot(Psix);

                        // N with RWG.div X RWG.div, kernelN = kernelA1
                        double LocalMatrixN_ii_jj = __ldg(&W[QudPt])                        //
                                                    * KernelA1(abc_alpha[4 * fieldIdx],     //
                                                               abc_alpha[4 * fieldIdx + 1], //
                                                               abc_alpha[4 * fieldIdx + 2], //
                                                               abc_alpha[4 * fieldIdx + 3], //
                                                               chi_tau, chi_t, chi_t - chi_tau) *
                                                    4 * fluxI * fluxJ;

                        // Combining all the atomic adds above into a single one
                        localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += 1 / (2 * mu0) * ((1 + mu / mu0) * (TnA[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj * TnA[DofsJ[permJ[jj]]]        // A1
                                                                                                                         + 2 * TnA[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj * TnA[DofsJ[permJ[jj]]]) // A2
                                                                                                       - 4 * (TnA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj * TdA[DofsJ[permJ[jj]]]                   // C1
                                                                                                              + TdA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj * TnA[DofsJ[permJ[jj]]]                 // C2
                                                                                                              + TnA[DofsI[permI[ii]]] * LocalMatrixC3_ii_jj * TdA[DofsJ[permJ[jj]]])                // C3
                                                                                                       + (1 + mu0 / mu) * -TdA[DofsI[permI[ii]]] * LocalMatrixN_ii_jj * TdA[DofsJ[permJ[jj]]])      // N
                                                                                      + mu / 2 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj * Mxn_coeffs[DofsJ[permJ[jj]]]                  // red_remaining
                                                                                      + mu * TnA[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj * Mxn_coeffs[DofsJ[permJ[jj]]]                             // red_l_M
                                                                                      + mu * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj * Mxn_coeffs[DofsJ[permJ[jj]]]                      // blue_remaining
                                                                                      + mu * (Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj * TnA[DofsJ[permJ[jj]]] +                          // blue_l_M
                                                                                              TnA[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj * Mxn_coeffs[DofsJ[permJ[jj]]]) +
                                                                                      -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj * TdA[DofsJ[permJ[jj]]]    // MC1
                                                                                      + -mu0 * TdA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj * Mxn_coeffs[DofsJ[permJ[jj]]]  // MC2
                                                                                      + -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC3_ii_jj * TdA[DofsJ[permJ[jj]]]; // MC3
                    }
                }
                // Combining all the atomic adds above into a single one
                // atomicAdd(&shapeDerivative[fieldIdx], localval);
            }
        }
    }
    // Writing to global memory
    for (int l = 0; l < Nabc_alpha; ++l)
    {
        atomicAdd(&shapeDerivative[l], localShapeDerivatives[l + threadIdx.x * Nabc_alpha]);
        // localShapeDerivatives[l + threadIdx.x * Nabc_alpha] = 0;
    }
}