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

__device__ Eigen::Vector3d Vel(Eigen::Vector3d X)
{
    // return Eigen::Vector3d(1, 0, 0);
    return Eigen::Vector3d(X(0) * X(1) * X(2), 0, 0);
}

__device__ Eigen::Matrix3d DVel(Eigen::Vector3d X)
{
    Eigen::Matrix3d out;
    out << X(1) * X(2), X(0) * X(2), X(0) * X(1), 0, 0, 0, 0, 0, 0;
    return out;
    // return Eigen::Matrix3d::Zero(3, 3);
}

__device__ double KernelA1(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    double znorm = YmX.norm();
    return YmX.dot(Vel(X) - Vel(Y)) / (4 * M_PI) / (znorm * znorm * znorm);
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

__device__ Eigen::Vector3d KernelC3(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    double znorm = YmX.norm();
    return -3. / (4 * M_PI) * YmX * YmX.dot(Vel(Y) - Vel(X)) / (znorm * znorm * znorm * znorm * znorm) + 1. / (4 * M_PI) * (Vel(Y) - Vel(X)) / (znorm * znorm * znorm);
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
                                       int NThreads, const int *I, const int *J, const int *relation,
                                       const double *W0, const double *X0, int Nq0,
                                       const double *W1, const double *X1, int Nq1,
                                       const double *W2, const double *X2, int Nq2,
                                       const double *W3, const double *X3, int Nq3,
                                       double *A1, double *A2, double *C1, double *C2, double *C3, double *N,
                                       const double *TdA, const double *TnA,
                                       const int *Elements, const double *Vertices, const double *Normals, const double *Areas,
                                       const int *Elt2DofTest, const int *Elt2DofTrial,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest)
/* double *testout, double *testABCi, double *testABCj,
int *orig_elti, int *orig_eltj, int *modif_elti, int *modif_eltj) */
{
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

    if (blockIdx.x == 0 && threadIdx.x == 0)
    {
        printf("Number of blocks: %d \n ", gridDim.x);
        printf("Threads per block: %d \n ", blockDim.x);
        printf("Total interactions: %d , Interactions per thread: %d \n ", NInteractions, InteractionsPerThread);
    }

    // Looping over all assigned interactions
    for (int idx = 0; idx < InteractionsPerThread; ++idx)
    {
        if (blockIdx.x == 0 && threadIdx.x == 0)
        {
            printf("In block 0 thread 0 computing interaction no. : %d \n ", idx);
        }
        // The interaction number
        // int InteractionIdx = ThreadID + NThreads * idx;
        int InteractionIdx = ThreadID * InteractionsPerThread + idx;

        // Preparing variables
        Eigen::Vector3d Ai, Bi, Ci, Aj, Bj, Cj;
        Eigen::MatrixXd Ei(3, 2), Ej(3, 2);

        int intersection[3], diffI[3], diffJ[3];

        const double *W = NULL;
        const double *X = NULL;
        int NQudPts = 0;

        // Make sure that the last thread stays in limit
        if (InteractionIdx >= NInteractions)
            break;

        // The pair of panels
        int i = I[InteractionIdx], j = J[InteractionIdx];

        /* printf("Interaction  (%d, %d) \n", i, j); */

        double g_tau = 2 * Areas[i], g_t = 2 * Areas[j];

        // Obtaining the normals
        /* Eigen::Vector3d normalx(Normals[i], Normals[i + NTriangles], Normals[i + 2 * NTriangles]);
        Eigen::Vector3d normaly(Normals[j], Normals[j + NTriangles], Normals[j + 2 * NTriangles]); */

        Eigen::Vector3d normalx(Normals[3 * i], Normals[3 * i + 1], Normals[3 * i + 2]);
        Eigen::Vector3d normaly(Normals[3 * j], Normals[3 * j + 1], Normals[3 * j + 2]);

        /* int EltI[] = {Elements[i],
                      Elements[i + NTriangles],
                      Elements[i + 2 * NTriangles]};

        int EltJ[] = {Elements[j],
                      Elements[j + NTriangles],
                      Elements[j + 2 * NTriangles]}; */

        int EltI[] = {Elements[3 * i],
                      Elements[3 * i + 1],
                      Elements[3 * i + 2]};

        int EltJ[] = {Elements[3 * j],
                      Elements[3 * j + 1],
                      Elements[3 * j + 2]};

        int origEltI[] = {EltI[0], EltI[1], EltI[2]};
        int origEltJ[] = {EltJ[0], EltJ[1], EltJ[2]};

        /* int DofsI[] = {Elt2DofTest[i], Elt2DofTest[i + NTriangles], Elt2DofTest[i + 2 * NTriangles]};
        int DofsJ[] = {Elt2DofTrial[j], Elt2DofTrial[j + NTriangles], Elt2DofTrial[j + 2 * NTriangles]}; */

        int DofsI[] = {Elt2DofTest[3 * i], Elt2DofTest[3 * i + 1], Elt2DofTest[3 * i + 2]};
        int DofsJ[] = {Elt2DofTrial[3 * j], Elt2DofTrial[3 * j + 1], Elt2DofTrial[3 * j + 2]};

        // Original permutation of elements
        int permI[] = {0, 1, 2};
        int permJ[] = {0, 1, 2};

        if (relation[InteractionIdx] == 0) // No interaction
        {
            // Computing Quadrature
            W = W0;
            X = X0;
            NQudPts = Nq0;
        }
        else if (relation[InteractionIdx] == 1) // Common vertex
        {
            IntersectionDiff(EltI, EltJ, intersection, diffI, diffJ);

            for (int l = 0; l < 3; ++l)
            {
                // Permutation for I
                if (EltI[l] == intersection[0])
                {
                    permI[0] = l;
                }
                else if (EltI[l] == diffI[0])
                {
                    permI[1] = l;
                }
                else if (EltI[l] == diffI[1])
                {
                    permI[2] = l;
                }

                // Permutation for J
                if (EltJ[l] == intersection[0])
                {
                    permJ[0] = l;
                }
                else if (EltJ[l] == diffJ[0])
                {
                    permJ[1] = l;
                }
                else if (EltJ[l] == diffJ[1])
                {
                    permJ[2] = l;
                }
            }

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

            for (int l = 0; l < 3; ++l)
            {
                // Permutation for I
                if (EltI[l] == intersection[0])
                {
                    permI[0] = l;
                }
                else if (EltI[l] == intersection[1])
                {
                    permI[1] = l;
                }
                else if (EltI[l] == diffI[0])
                {
                    permI[2] = l;
                }

                // Permutation for J
                if (EltJ[l] == intersection[0])
                {
                    permJ[0] = l;
                }
                else if (EltJ[l] == intersection[1])
                {
                    permJ[1] = l;
                }
                else if (EltJ[l] == diffJ[0])
                {
                    permJ[2] = l;
                }
            }

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
        }

        // Vertices of element i
        Ai = Eigen::Vector3d(Vertices[3 * EltI[0]], Vertices[3 * EltI[0] + 1], Vertices[3 * EltI[0] + 2]);
        Bi = Eigen::Vector3d(Vertices[3 * EltI[1]], Vertices[3 * EltI[1] + 1], Vertices[3 * EltI[1] + 2]);
        Ci = Eigen::Vector3d(Vertices[3 * EltI[2]], Vertices[3 * EltI[2] + 1], Vertices[3 * EltI[2] + 2]);

        // Vertices of element j
        Aj = Eigen::Vector3d(Vertices[3 * EltJ[0]], Vertices[3 * EltJ[0] + 1], Vertices[3 * EltJ[0] + 2]);
        Bj = Eigen::Vector3d(Vertices[3 * EltJ[1]], Vertices[3 * EltJ[1] + 1], Vertices[3 * EltJ[1] + 2]);
        Cj = Eigen::Vector3d(Vertices[3 * EltJ[2]], Vertices[3 * EltJ[2] + 1], Vertices[3 * EltJ[2] + 2]);

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

        Eigen::Matrix3d LocalMatrix = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrixA1 = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrixA2 = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrixC1 = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrixC3 = Eigen::MatrixX3d::Zero(3, 3);
        Eigen::Matrix3d LocalMatrixN = Eigen::MatrixX3d::Zero(3, 3);

        // P0 X P0
        /* for (int ii = 0; ii < NRSFTest; ++ii)
        {
            double Psix = g_tau;
            for (int jj = 0; jj < NRSFTrial; ++jj)
            {
                double Psiy = g_t;
                for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                {
                    Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[QudPt] + Ei.col(1) * X[QudPt + NQudPts];
                    Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[QudPt + 2 * NQudPts] + Ej.col(1) * X[QudPt + 3 * NQudPts];
                    LocalMatrix(ii, jj) += W[QudPt] * Psix * SLKernel(chi_tau, chi_t, chi_t - chi_tau) * Psiy;
                }
                GalerkinMatrix[i + TestDim * j] += LocalMatrix(ii, jj);

                // Atomic update of the galerkin matrix
                // double contribution = LocalMatrix(ii, jj);
                // atomicAdd(&GalerkinMatrix[i + TestDim * j], contribution);
            }
        } */

        // P1 X P1
        /* for (int ii = 0; ii < NRSFTest; ++ii)
        {

            for (int jj = 0; jj < NRSFTrial; ++jj)
            {

                for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                {
                    Eigen::Vector3d RSFsX(1 - X[4 * QudPt] - X[4 * QudPt + 1], X[4 * QudPt], X[4 * QudPt + 1]);
                    Eigen::Vector3d RSFsY(1 - X[4 * QudPt + 2] - X[4 * QudPt + 3], X[4 * QudPt + 2], X[4 * QudPt + 3]);

                    RSFsX *= g_tau;
                    RSFsY *= g_t;

                    Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
                    Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[4 * QudPt + 2] + Ej.col(1) * X[4 * QudPt + 3];

                    LocalMatrix(ii, jj) += W[QudPt] * RSFsX(ii) * SLKernel(chi_tau, chi_t, chi_t - chi_tau) * RSFsY(jj);
                }
                // GalerkinMatrix[EltI[ii] + TestDim * EltJ[jj]] += LocalMatrix(ii, jj);

                //  Atomic update of the galerkin matrix
                atomicAdd(&GalerkinMatrix[EltI[ii] + TestDim * EltJ[jj]], LocalMatrix(ii, jj));
            }
        } */

        // P0 X ntimes(P1)
        /* for (int ii = 0; ii < NRSFTest; ++ii)
        {
            double Psix = g_tau;

            for (int jj = 0; jj < NRSFTrial; ++jj)
            {

                for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                {

                    Eigen::Vector3d RSFsY(1 - X[QudPt + 2 * NQudPts] - X[QudPt + 3 * NQudPts], X[QudPt + 2 * NQudPts], X[QudPt + 3 * NQudPts]);

                    RSFsY *= g_t;

                    Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[QudPt] + Ei.col(1) * X[QudPt + NQudPts];
                    Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[QudPt + 2 * NQudPts] + Ej.col(1) * X[QudPt + 3 * NQudPts];

                    LocalMatrix(ii, jj) += W[QudPt] * Psix * DLKernel(chi_tau, chi_t, chi_t - chi_tau).dot(normaly) * RSFsY(jj);
                }
                // GalerkinMatrix[i + TestDim * EltJ[jj]] += LocalMatrix(ii, jj);

                //  Atomic update of the galerkin matrix
                atomicAdd(&GalerkinMatrix[i + TestDim * EltJ[jj]], LocalMatrix(ii, jj));
            }
        } */

        // RWG X RWG
        for (int ii = 0; ii < 3; ++ii)
        {
            int iip1 = (permI[ii] + 1) % 3;
            int iip2 = (iip1 + 1) % 3;

            double fluxI = origEltI[iip1] < origEltI[iip2] ? 1. : -1.;
            // double RWGX_ref_0 = -ii % 2, RWGX_ref_1 = -ii / 2;

            for (int jj = 0; jj < 3; ++jj)
            {
                int jjp1 = (permJ[jj] + 1) % 3;
                int jjp2 = (jjp1 + 1) % 3;

                double fluxJ = origEltJ[jjp1] < origEltJ[jjp2] ? 1. : -1.;
                // double RWGY_ref_0 = -jj % 2, RWGY_ref_1 = -jj / 2;

                for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                {
                    /* if (blockIdx.x == 0 && threadIdx.x == 0)
                    {
                        printf("Qud pt %d\n", QudPt);
                    } */
                    // Reference basis RT0
                    /* Eigen::MatrixXd RWGX_ref(3, 2); // Rows represent the 3 RSFs
                    RWGX_ref << X[4 * QudPt], X[4 * QudPt + 1],
                        X[4 * QudPt] - 1, X[4 * QudPt + 1],
                        X[4 * QudPt], X[4 * QudPt + 1] - 1;

                    Eigen::MatrixXd RWGY_ref(3, 2); // Rows represent the 3 RSFs
                    RWGY_ref << X[4 * QudPt + 2], X[4 * QudPt + 3],
                        X[4 * QudPt + 2] - 1, X[4 * QudPt + 3],
                        X[4 * QudPt + 2], X[4 * QudPt + 3] - 1;

                    Eigen::Vector3d Psix = fluxI * Ei * RWGX_ref.row(ii).transpose();
                    Eigen::Vector3d Psiy = fluxJ * Ej * RWGY_ref.row(jj).transpose(); */

                    double RWGX_ref_0 = X[4 * QudPt] - ii % 2;
                    double RWGX_ref_1 = X[4 * QudPt + 1] - ii / 2;

                    double RWGY_ref_0 = X[4 * QudPt + 2] - jj % 2;
                    double RWGY_ref_1 = X[4 * QudPt + 3] - jj / 2;

                    // RWG elements
                    Eigen::Vector3d Psix = fluxI * (Ei.col(0) * RWGX_ref_0 + Ei.col(1) * RWGX_ref_1);
                    Eigen::Vector3d Psiy = fluxJ * (Ej.col(0) * RWGY_ref_0 + Ej.col(1) * RWGY_ref_1);

                    Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
                    Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[4 * QudPt + 2] + Ej.col(1) * X[4 * QudPt + 3];

                    // A1 with RWG . RWG
                    LocalMatrixA1(ii, jj) += W[QudPt] * KernelA1(chi_tau, chi_t, chi_t - chi_tau) * Psiy.dot(Psix);

                    // A2 with DVelRWG(y) . RWG(x)
                    LocalMatrixA2(ii, jj) += W[QudPt] * KernelA2(chi_tau, chi_t, chi_t - chi_tau) * (DVel(chi_t) * Psiy).dot(Psix);

                    // C1 with DVelRWG(y) X RWG(X)
                    LocalMatrixC1(ii, jj) += W[QudPt] * (KernelC1(chi_tau, chi_t, chi_t - chi_tau).cross(DVel(chi_t) * Psiy)).dot(Psix);

                    // C3 with RWG X RWG
                    LocalMatrixC3(ii, jj) += W[QudPt] * (KernelC3(chi_tau, chi_t, chi_t - chi_tau).cross(Psiy)).dot(Psix);

                    // N with RWG.div X RWG.div, kernelN = kernelA1
                    LocalMatrixN(ii, jj) += W[QudPt] * KernelA1(chi_tau, chi_t, chi_t - chi_tau) * 4 * fluxI * fluxJ;
                }
                // GalerkinMatrix[i + TestDim * EltJ[jj]] += LocalMatrix(ii, jj);

                // Accumulating the shape derivative
                // Local matrix (ii,jj) contains part of the global element DofsI[permI[ii]] , DofsJ[permJ[jj]]

                // TnA' * A1mat * TnA
                atomicAdd(A1, TnA[DofsI[permI[ii]]] * LocalMatrixA1(ii, jj) * TnA[DofsJ[permJ[jj]]]);

                // 2 * TnA' * A2mat * TnA
                atomicAdd(A2, 2 * TnA[DofsI[permI[ii]]] * LocalMatrixA2(ii, jj) * TnA[DofsJ[permJ[jj]]]);

                // C1 = TnA' * C1mat * TdA
                atomicAdd(C1, TnA[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * TdA[DofsJ[permJ[jj]]]);

                // C2 = TdA' * C1mat * TnA
                atomicAdd(C2, TdA[DofsI[permI[ii]]] * LocalMatrixC1(ii, jj) * TnA[DofsJ[permJ[jj]]]);

                // C3 = TnA' * C3mat * TdA
                atomicAdd(C3, TnA[DofsI[permI[ii]]] * LocalMatrixC3(ii, jj) * TdA[DofsJ[permJ[jj]]]);

                // N = -TdA' * Nmat * TdA
                atomicAdd(N, -TdA[DofsI[permI[ii]]] * LocalMatrixN(ii, jj) * TdA[DofsJ[permJ[jj]]]);

                //  Atomic update of the galerkin matrix
                // atomicAdd(&GalerkinMatrix[DofsI[permI[ii]] + TestDim * DofsJ[permJ[jj]]], LocalMatrix(ii, jj));
            }
        }
    }
}