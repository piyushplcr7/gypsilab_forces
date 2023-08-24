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

__device__ double KV(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
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

__device__ double KernelOld(int a, int b, int c, int alpha, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    // 1/4pi * z.(vel(x) - vel(y))/||z||^3
    return YmX.dot(Vel(a, b, c, alpha, X) - Vel(a, b, c, alpha, Y)) / (4 * M_PI) / (znorm * znorm * znorm);
}

__device__ Eigen::Vector3d KernelIntegrable(int a, int b, int c, int alpha, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    // 3/4pi * z.(vel(y) - vel(x))/||z||^5  z (.n handled outside)
    return 3 * YmX.dot(Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) / (4 * M_PI) / pow(znorm, 5) * YmX;
}

__device__ Eigen::Vector3d KernelComb(int a, int b, int c, int alpha, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();

    return 1. / (4 * M_PI) * (-DVel(a, b, c, alpha, Y) * YmX + Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) / pow(znorm, 3);
}

__global__ void computeShapeDerivative(int TrialDim, int TestDim, int NTriangles, int NVertices, int NInteractions,
                                       int NThreads, const unsigned short int *I, const unsigned short int *J, const int *relation,
                                       const double *W0, const double *X0, int Nq0,
                                       const double *W1, const double *X1, int Nq1,
                                       const double *W2, const double *X2, int Nq2,
                                       const double *W3, const double *X3, int Nq3,
                                       double *t7, double *t1nt2, double *t56nt6, double *shapeDerivative,
                                       const double *Tdu, const double *HJn_coeffs,
                                       const int *Elements, const double *Vertices, const double *Normals, const double *Areas,
                                       const int *Elt2DofTest, const int *Elt2DofTrial,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest,
                                       const int *abc_alpha, int Nabc_alpha,
                                       const int *permII, const int *permJJ)
{
    __shared__ double localShapeDerivatives[81 * 32];
    // Initialization of shared memory
    for (int l = 0; l < Nabc_alpha; ++l)
    {
        localShapeDerivatives[l + threadIdx.x * Nabc_alpha] = 0;
    }

    int ThreadID = blockIdx.x * blockDim.x + threadIdx.x;

    // Matrix size is NTriangles^2. Each entry is assigned to a thread
    // InteractionsPerThread gives the max no.of element interactions assigned
    // to a thread
    int InteractionsPerThread = ceil(double(NInteractions) / double(NThreads));

    // Looping over all assigned interactions
    for (int idx = 0; idx < InteractionsPerThread; ++idx)
    {
        // The interaction number
        int InteractionIdx = ThreadID * InteractionsPerThread + idx;

        // Preparing variables
        Eigen::Vector3d Ai, Bi, Ci, Aj, Bj, Cj;
        Eigen::Vector3d Eicol0, Eicol1, Ejcol0, Ejcol1;

        // Make sure that the last thread stays in limit
        if (InteractionIdx >= NInteractions)
            break;

        // The pair of panels
        int i = __ldg(&I[InteractionIdx]), j = __ldg(&J[InteractionIdx]);
        double g_tau = 2 * __ldg(&Areas[i]), g_t = 2 * __ldg(&Areas[j]);

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
        Eicol0 = Bi - Ai;
        Eicol1 = Ci - Ai;

        Ejcol0 = Bj - Aj;
        Ejcol1 = Cj - Aj;

        Eigen::Matrix2d EtEi;
        // EtEi = Ei.transpose() * Ei;
        EtEi << Eicol0.dot(Eicol0), Eicol0.dot(Eicol1), Eicol1.dot(Eicol0), Eicol1.dot(Eicol1);
        Eigen::Matrix2d EtEj;
        // EtEj = Ej.transpose() * Ej;
        EtEj << Ejcol0.dot(Ejcol0), Ejcol0.dot(Ejcol1), Ejcol1.dot(Ejcol0), Ejcol1.dot(Ejcol1);

        double deti = EtEi(0, 0) * EtEi(1, 1) - EtEi(0, 1) * EtEi(1, 0);
        double detj = EtEj(0, 0) * EtEj(1, 1) - EtEj(0, 1) * EtEj(1, 0);

        // Dxyi and Dxyj are the inverses of EtEi and EtEj respectively
        Eigen::Matrix2d Dxyi = -EtEi;
        Eigen::Matrix2d Dxyj = -EtEj;

        Dxyi(0, 0) = EtEi(1, 1);
        Dxyi(1, 1) = EtEi(0, 0);

        Dxyj(0, 0) = EtEj(1, 1);
        Dxyj(1, 1) = EtEj(0, 0);

        Dxyi /= deti;
        Dxyj /= detj;

        // Eigen::MatrixXd DCVi = Ei * Dxyi, DCVj = Ej * Dxyj;
        //  Manually computing the above products
        Eigen::Vector3d DCVicol0 = Eicol0 * Dxyi(0, 0) + Eicol1 * Dxyi(0, 1);
        Eigen::Vector3d DCVicol1 = Eicol0 * Dxyi(1, 0) + Eicol1 * Dxyi(1, 1);

        Eigen::Vector3d DCVjcol0 = Ejcol0 * Dxyj(0, 0) + Ejcol1 * Dxyj(0, 1);
        Eigen::Vector3d DCVjcol1 = Ejcol0 * Dxyj(1, 0) + Ejcol1 * Dxyj(1, 1);

        // Mixed P0XP0, P1XP1, P0XP1
        for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
        {
            // P0 X P0
            double Psix_P0 = g_tau;
            double Psiy_P0 = g_t;

            Eigen::Vector3d chi_tau = Ai + Eicol0 * __ldg(&X[4 * QudPt]) + Eicol1 * __ldg(&X[4 * QudPt + 1]);
            Eigen::Vector3d chi_t = Aj + Ejcol0 * __ldg(&X[4 * QudPt + 2]) + Ejcol1 * __ldg(&X[4 * QudPt + 3]);

            for (int fieldIdx = 0; fieldIdx < Nabc_alpha; ++fieldIdx)
            {
                double Local_kerneloldmat_P0_P0 = __ldg(&W[QudPt]) * Psix_P0 * KernelOld(__ldg(&abc_alpha[4 * fieldIdx]),     //
                                                                                         __ldg(&abc_alpha[4 * fieldIdx + 1]), //
                                                                                         __ldg(&abc_alpha[4 * fieldIdx + 2]), //
                                                                                         __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau) *
                                                  Psiy_P0;

                localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += -0.5 * HJn_coeffs[i] * Local_kerneloldmat_P0_P0 * HJn_coeffs[j];

                // Double loop for P1
                // Trial
                for (int jj = 0; jj < 3; ++jj)
                {
                    // no. jj rsf
                    int Psiy_nxgradP1_0 = (3 - jj) % 3 - 1;
                    int Psiy_nxgradP1_1 = jj - 1;

                    Eigen::Vector3d Psiy_nxgradP1 = g_t * normaly.cross(DCVjcol0 * Psiy_nxgradP1_0 + DCVjcol1 * Psiy_nxgradP1_1);

                    Eigen::Vector3d RSFsY_P1(1 - __ldg(&X[4 * QudPt + 2]) - __ldg(&X[4 * QudPt + 3]), __ldg(&X[4 * QudPt + 2]), __ldg(&X[4 * QudPt + 3]));
                    RSFsY_P1 *= g_t;

                    // P0 X ntimes(P1) (test X trial)
                    double LocalMatrix_kernelintegrablemat_jj = __ldg(&W[QudPt]) * KernelIntegrable(__ldg(&abc_alpha[4 * fieldIdx]),     //
                                                                                                    __ldg(&abc_alpha[4 * fieldIdx + 1]), //
                                                                                                    __ldg(&abc_alpha[4 * fieldIdx + 2]), //
                                                                                                    __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau)
                                                                                       .dot(normaly) *
                                                                RSFsY_P1(jj) * Psix_P0;

                    double LocalMatrix_combkernelmat_jj = __ldg(&W[QudPt]) * KernelComb(__ldg(&abc_alpha[4 * fieldIdx]),     //
                                                                                        __ldg(&abc_alpha[4 * fieldIdx + 1]), //
                                                                                        __ldg(&abc_alpha[4 * fieldIdx + 2]), //
                                                                                        __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau)
                                                                                 .dot(normaly) *
                                                          RSFsY_P1(jj) * Psix_P0;

                    localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += HJn_coeffs[i] * (-LocalMatrix_kernelintegrablemat_jj + LocalMatrix_combkernelmat_jj) * Tdu[EltJ[jj]];

                    // Test
                    for (int ii = 0; ii < 3; ++ii)
                    {
                        // no. jj rsf
                        int Psix_nxgradP1_0 = (3 - ii) % 3 - 1;
                        int Psix_nxgradP1_1 = ii - 1;

                        Eigen::Vector3d Psix_nxgradP1 = g_tau * normalx.cross(DCVicol0 * Psix_nxgradP1_0 + DCVicol1 * Psix_nxgradP1_1);

                        Eigen::Vector3d RSFsX_P1(1 - __ldg(&X[4 * QudPt]) - __ldg(&X[4 * QudPt + 1]), __ldg(&X[4 * QudPt]), __ldg(&X[4 * QudPt + 1]));
                        RSFsX_P1 *= g_tau;

                        // KernelOld with nxgradP1 X nxgradP1
                        double LocalMatrix_kerneloldmat_nxgradP1_nxgradP1_ii_jj = __ldg(&W[QudPt]) * KernelOld(__ldg(&abc_alpha[4 * fieldIdx]),     //
                                                                                                               __ldg(&abc_alpha[4 * fieldIdx + 1]), //
                                                                                                               __ldg(&abc_alpha[4 * fieldIdx + 2]), //
                                                                                                               __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_tau, chi_t, chi_t - chi_tau) *
                                                                                  Psiy_nxgradP1.dot(Psix_nxgradP1);
                        double LocalMatrix_SL_Dvelnxgrad_nxgrad_ii_jj = __ldg(&W[QudPt]) * KV(chi_tau, chi_t, chi_t - chi_tau) * (DVel(__ldg(&abc_alpha[4 * fieldIdx]),     //
                                                                                                                                       __ldg(&abc_alpha[4 * fieldIdx + 1]), //
                                                                                                                                       __ldg(&abc_alpha[4 * fieldIdx + 2]), //
                                                                                                                                       __ldg(&abc_alpha[4 * fieldIdx + 3]), chi_t) *
                                                                                                                                  Psiy_nxgradP1)
                                                                                                                                     .dot(Psix_nxgradP1);

                        localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += Tdu[EltI[ii]] * (0.5 * LocalMatrix_kerneloldmat_nxgradP1_nxgradP1_ii_jj + LocalMatrix_SL_Dvelnxgrad_nxgrad_ii_jj) * Tdu[EltJ[jj]];
                    }
                }
            }
        }
    }

    // Writing to global memory
    for (int l = 0; l < Nabc_alpha; ++l)
    {
        atomicAdd(&shapeDerivative[l], localShapeDerivatives[l + threadIdx.x * Nabc_alpha]);
    }
}