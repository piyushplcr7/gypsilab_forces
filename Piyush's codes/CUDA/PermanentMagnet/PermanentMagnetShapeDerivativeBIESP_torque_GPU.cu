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

__device__ Eigen::Vector3d Vel(int dir, Eigen::Vector3d Xcg, Eigen::Vector3d X)
{
    Eigen::Vector3d axis = Eigen::Vector3d::Zero(3);
    axis(dir) = 1;
    Eigen::Vector3d out = axis.cross(X - Xcg);
    return out;
}

__device__ Eigen::Matrix3d DVel(int dir, Eigen::Vector3d Xcg, Eigen::Vector3d X)
{
    Eigen::Vector3d axis = Eigen::Vector3d::Zero(3);
    axis(dir) = 1;
    Eigen::Matrix3d out(3, 3);
    out << 0, -axis(2), axis(1),
        axis(2), 0, -axis(0),
        -axis(1), axis(0), 0;
    return out;
}

__device__ double KernelOld(int dir, Eigen::Vector3d Xcg, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    // 1/4pi * z.(vel(x) - vel(y))/||z||^3
    return YmX.dot(Vel(dir, Xcg, X) - Vel(dir, Xcg, Y)) / (4 * M_PI) / (znorm * znorm * znorm);
}

__device__ Eigen::Vector3d KernelIntegrable(int dir, Eigen::Vector3d Xcg, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    // 3/4pi * z.(vel(y) - vel(x))/||z||^5  z (.n handled outside)
    return 3 * YmX.dot(Vel(dir, Xcg, Y) - Vel(dir, Xcg, X)) / (4 * M_PI) / pow(znorm, 5) * YmX;
}

__device__ Eigen::Vector3d KernelComb(int dir, Eigen::Vector3d Xcg, Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();

    return 1. / (4 * M_PI) * (-DVel(dir, Xcg, Y) * YmX + Vel(dir, Xcg, Y) - Vel(dir, Xcg, X)) / pow(znorm, 3);
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
                                       double *dbv_ds, double *dbw_ds, double *dbk_ds, double *l1, double *l2, double *r1,
                                       double *shapeDerivative,
                                       const double *Tdu, const double *Tnu, const double *lambda_coeffs, double mu0,
                                       const int *Elements, const double *Vertices, const double *Normals, const double *Areas,
                                       const int *Elt2DofTest, const int *Elt2DofTrial,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest,
                                       const int *permII, const int *permJJ,
                                       const double *Xcg)
/* double *testout, double *testABCi, double *testABCj,
int *orig_elti, int *orig_eltj, int *modif_elti, int *modif_eltj) */
{
    // shared memory for writing shape derivatives
    // Need to specify the size beforehand, hard coding Nabc_alpha = 81, blockDim.x = 32
    __shared__ double localShapeDerivatives[3 * 32];
    int Nabc_alpha = 3;
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
        // Eigen::MatrixXd Ei(3, 2), Ej(3, 2);
        Eigen::Vector3d Eicol0, Eicol1, Ejcol0, Ejcol1;

        int intersection[3], diffI[3], diffJ[3];

        /* const double *W = NULL;
        const double *X = NULL;
        int NQudPts = 0; */

        // Make sure that the last thread stays in limit
        if (InteractionIdx >= NInteractions)
            break;

        // The pair of panels
        int i = __ldg(&I[InteractionIdx]), j = __ldg(&J[InteractionIdx]);

        /* printf("Interaction  (%d, %d) \n", i, j); */

        double g_tau = 2 * __ldg(&Areas[i]), g_t = 2 * __ldg(&Areas[j]);

        // Obtaining the normals
        /* Eigen::Vector3d normalx(Normals[i], Normals[i + NTriangles], Normals[i + 2 * NTriangles]);
        Eigen::Vector3d normaly(Normals[j], Normals[j + NTriangles], Normals[j + 2 * NTriangles]); */

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

        /* if (relation[InteractionIdx] == 0) // No interaction
        {
            // Computing Quadrature
            W = W0;
            X = X0;
            NQudPts = Nq0;
        }
        else if (relation[InteractionIdx] == 1) // Common vertex
        {
            IntersectionDiff(EltI, EltJ, intersection, diffI, diffJ);

            // Changing EltI into ABCI
            EltI[0] = intersection[0];
            EltI[1] = diffI[0];
            EltI[2] = diffI[1];

            // Changing EltI into ABCJ
            EltJ[0] = intersection[0];
            EltJ[1] = diffJ[0];
            EltJ[2] = diffJ[1];

            // Computing Quadrature
            W = W1;
            X = X1;
            NQudPts = Nq1;
        }
        else if (relation[InteractionIdx] == 2) // Common edge
        {
            IntersectionDiff(EltI, EltJ, intersection, diffI, diffJ);

            EltI[0] = intersection[0];
            EltI[1] = intersection[1];
            EltI[2] = diffI[0];

            EltJ[0] = intersection[0];
            EltJ[1] = intersection[1];
            EltJ[2] = diffJ[0];

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

        /* Ei.col(0) = Bi - Ai;
        Ei.col(1) = Ci - Ai;

        Ej.col(0) = Bj - Aj;
        Ej.col(1) = Cj - Aj; */

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

        /* Eigen::Matrix3d LocalMatrix_kerneloldmat_nxgradP1_nxgradP1 = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrix_SL_Dvelnxgrad_nxgrad = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Vector3d LocalMatrix_kernelintegrablemat = Eigen::Vector3d::Zero(3);
        Eigen::Vector3d LocalMatrix_combkernelmat = Eigen::Vector3d::Zero(3);
        double Local_kerneloldmat_P0_P0 = 0; */

        Eigen::Vector3d XCG;
        XCG(0) = Xcg[0];
        XCG(1) = Xcg[1];
        XCG(2) = Xcg[2];

        // Mixed P0XP0, P1XP1, P0XP1
        for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
        {
            // P0 X P0
            double Psix_P0 = g_tau;
            double Psiy_P0 = g_t;

            Eigen::Vector3d chi_tau = Ai + Eicol0 * __ldg(&X[4 * QudPt]) + Eicol1 * __ldg(&X[4 * QudPt + 1]);
            Eigen::Vector3d chi_t = Aj + Ejcol0 * __ldg(&X[4 * QudPt + 2]) + Ejcol1 * __ldg(&X[4 * QudPt + 3]);

            // Looping over different fields
            for (int fieldIdx = 0; fieldIdx < Nabc_alpha; ++fieldIdx)
            {
                double Local_kerneloldmat_P0_P0 = __ldg(&W[QudPt]) * Psix_P0 * KernelOld(fieldIdx, XCG, //
                                                                                         chi_tau, chi_t, chi_t - chi_tau) *
                                                  Psiy_P0;

                localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += mu0 * (-Tnu[i] * Local_kerneloldmat_P0_P0 * Tnu[j])                          // dbv_ds
                                                                              + mu0 * Tnu[i] * Local_kerneloldmat_P0_P0 * lambda_coeffs[j]                 // l1
                                                                              + -mu0 / 2 * lambda_coeffs[i] * Local_kerneloldmat_P0_P0 * lambda_coeffs[j]; // r1

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
                    double LocalMatrix_kernelintegrablemat_jj = __ldg(&W[QudPt]) * KernelIntegrable(fieldIdx, XCG, //
                                                                                                    chi_tau, chi_t, chi_t - chi_tau)
                                                                                       .dot(normaly) *
                                                                RSFsY_P1(jj) * Psix_P0;

                    double LocalMatrix_combkernelmat_jj = __ldg(&W[QudPt]) * KernelComb(fieldIdx, XCG, //
                                                                                        chi_tau, chi_t, chi_t - chi_tau)
                                                                                 .dot(normaly) *
                                                          RSFsY_P1(jj) * Psix_P0;

                    localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += 2 * mu0 * Tnu[i] * (LocalMatrix_kernelintegrablemat_jj - LocalMatrix_combkernelmat_jj) * Tdu[EltJ[jj]]          // dbk_ds
                                                                                  - mu0 * lambda_coeffs[i] * (LocalMatrix_kernelintegrablemat_jj - LocalMatrix_combkernelmat_jj) * Tdu[EltJ[jj]]; // l2

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
                        double LocalMatrix_kerneloldmat_nxgradP1_nxgradP1_ii_jj = __ldg(&W[QudPt]) * KernelOld(fieldIdx, XCG, //
                                                                                                               chi_tau, chi_t, chi_t - chi_tau) *
                                                                                  Psiy_nxgradP1.dot(Psix_nxgradP1);
                        double LocalMatrix_SL_Dvelnxgrad_nxgrad_ii_jj = __ldg(&W[QudPt]) * KV(chi_tau, chi_t, chi_t - chi_tau) //
                                                                        * (DVel(fieldIdx, XCG, chi_t) *
                                                                           Psiy_nxgradP1)
                                                                              .dot(Psix_nxgradP1);

                        localShapeDerivatives[fieldIdx + threadIdx.x * Nabc_alpha] += mu0 * Tdu[EltI[ii]] * (LocalMatrix_kerneloldmat_nxgradP1_nxgradP1_ii_jj // dbw_ds
                                                                                                             + 2 * LocalMatrix_SL_Dvelnxgrad_nxgrad_ii_jj) *
                                                                                      Tdu[EltJ[jj]];
                    }
                }
            }
        }

        /* double shapeDerivativeForThisInteraction = 0;
        // Updating the outputs atomically
        shapeDerivativeForThisInteraction += mu0 * (-Tnu[i] * Local_kerneloldmat_P0_P0 * Tnu[j])                          // dbv_ds
                                             + mu0 * Tnu[i] * Local_kerneloldmat_P0_P0 * lambda_coeffs[j]                 // l1
                                             + -mu0 / 2 * lambda_coeffs[i] * Local_kerneloldmat_P0_P0 * lambda_coeffs[j]; // r1

        for (int jj = 0; jj < 3; ++jj)
        {
            shapeDerivativeForThisInteraction += 2 * mu0 * Tnu[i] * (LocalMatrix_kernelintegrablemat(jj) - LocalMatrix_combkernelmat(jj)) * Tdu[EltJ[jj]]          // dbk_ds
                                                 - mu0 * lambda_coeffs[i] * (LocalMatrix_kernelintegrablemat(jj) - LocalMatrix_combkernelmat(jj)) * Tdu[EltJ[jj]]; // l2

            for (int ii = 0; ii < 3; ++ii)
            {
                shapeDerivativeForThisInteraction += mu0 * Tdu[EltI[ii]] * (LocalMatrix_kerneloldmat_nxgradP1_nxgradP1(ii, jj) // dbw_ds
                                                                            + 2 * LocalMatrix_SL_Dvelnxgrad_nxgrad(ii, jj)) *
                                                     Tdu[EltJ[jj]];
            }
        } */

        /* atomicAdd(shapeDerivative, shapeDerivativeForThisInteraction); */
    }

    // Writing to global memory
    for (int l = 0; l < Nabc_alpha; ++l)
    {
        atomicAdd(&shapeDerivative[l], localShapeDerivatives[l + threadIdx.x * Nabc_alpha]);
        // localShapeDerivatives[l + threadIdx.x * Nabc_alpha] = 0;
    }
}