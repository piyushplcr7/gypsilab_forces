#include <cuda_runtime.h>
#include <stdio.h>
#include <cuda.h>
#include <math.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

__device__ double Kernel(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    // return 1;
    return 1. / (4 * M_PI) / YmX.norm();
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

// Defining the kernel for SS computation
/*
 * inputs
 *
 * NTriangles : Number of elements
 */
__global__ void computeShapeDerivative(int TrialDim, int TestDim, int NTriangles, int NVertices,
                                       int NThreads, const int *I, const int *J, const int *relation,
                                       const double *W0, const double *X0, int Nq0,
                                       const double *W1, const double *X1, int Nq1,
                                       const double *W2, const double *X2, int Nq2,
                                       const double *W3, const double *X3, int Nq3,
                                       double *shapeDerivative,
                                       double *GalerkinMatrix,
                                       const double *trial_vec, const double *test_vec,
                                       const int *Elements, const double *Vertices, const double *Normals, const double *Areas,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest)
/* double *testout, double *testABCi, double *testABCj,
int *orig_elti, int *orig_eltj, int *modif_elti, int *modif_eltj) */
{
    int ThreadID = blockIdx.x * blockDim.x + threadIdx.x;

    *shapeDerivative = 3.145;

    /* int idebug = 1, jdebug = 2; */

    // Size of the matrix
    // int NInteractions = TrialDim * TestDim;

    // Number of element interaction
    int NInteractions = NTriangles * NTriangles;

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

    // Stores the contribution to the shape derivative from the
    // local matrix that will be computed in this thread
    // double localval = 0;

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
        Eigen::Vector3d Ai, Bi, Ci, ABCi, permi, Aj, Bj, Cj, ABCj, permj;
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

        int EltI[] = {Elements[i],
                      Elements[i + NTriangles],
                      Elements[i + 2 * NTriangles]};

        int EltJ[] = {Elements[j],
                      Elements[j + NTriangles],
                      Elements[j + 2 * NTriangles]};

        // Storing the original elements
        /* if (i == idebug && j == jdebug)
        {
            // Storing the original element assigned at this point
            for (int idx = 0; idx < 3; ++idx)
            {
                orig_elti[idx] = EltI[idx];
                orig_eltj[idx] = EltJ[idx];
            }
        } */

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

            EltI[0] = intersection[0];
            EltI[1] = diffI[0];
            EltI[2] = diffI[1];

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
        }

        // Vertices of element i
        Ai = Eigen::Vector3d(Vertices[EltI[0]], Vertices[EltI[0] + NVertices], Vertices[EltI[0] + 2 * NVertices]);
        Bi = Eigen::Vector3d(Vertices[EltI[1]], Vertices[EltI[1] + NVertices], Vertices[EltI[1] + 2 * NVertices]);
        Ci = Eigen::Vector3d(Vertices[EltI[2]], Vertices[EltI[2] + NVertices], Vertices[EltI[2] + 2 * NVertices]);

        // Vertices of element j
        Aj = Eigen::Vector3d(Vertices[EltJ[0]], Vertices[EltJ[0] + NVertices], Vertices[EltJ[0] + 2 * NVertices]);
        Bj = Eigen::Vector3d(Vertices[EltJ[1]], Vertices[EltJ[1] + NVertices], Vertices[EltJ[1] + 2 * NVertices]);
        Cj = Eigen::Vector3d(Vertices[EltJ[2]], Vertices[EltJ[2] + NVertices], Vertices[EltJ[2] + 2 * NVertices]);

        // Storing ABC after reassignment from intersection and diff
        /* if (i == idebug && j == jdebug)
        {
            // storing the vertices using column major format

            for (int idx = 0; idx < 3; ++idx)
            {
                testABCi[0 + 3 * idx] = Ai(idx);
                testABCi[1 + 3 * idx] = Bi(idx);
                testABCi[2 + 3 * idx] = Ci(idx);

                testABCj[0 + 3 * idx] = Aj(idx);
                testABCj[1 + 3 * idx] = Bj(idx);
                testABCj[2 + 3 * idx] = Cj(idx);
            }
        }

        // Storing the modified element
        if (i == idebug && j == jdebug)
        {
            for (int idx = 0; idx < 3; ++idx)
            {
                modif_elti[idx] = EltI[idx];
                modif_eltj[idx] = EltJ[idx];
            }
        } */

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

        Eigen::MatrixXd LocalMatrix = Eigen::MatrixXd::Zero(NRSFTest, NRSFTrial);
        // Performing integration according to the spaces
        // P0
        if (TrialSpace == 0)
        {
            for (int ii = 0; ii < NRSFTest; ++ii)
            {
                double Psix = g_tau;
                for (int jj = 0; jj < NRSFTrial; ++jj)
                {
                    double Psiy = g_t;
                    for (int QudPt = 0; QudPt < NQudPts; ++QudPt)
                    {
                        Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[QudPt] + Ei.col(1) * X[QudPt + NQudPts];
                        Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[QudPt + 2 * NQudPts] + Ej.col(1) * X[QudPt + 3 * NQudPts];
                        LocalMatrix(ii, jj) += W[QudPt] * Psix * Kernel(chi_tau, chi_t, chi_t - chi_tau) * Psiy;

                        // Storing Kernel values
                        /* if (i == idebug && j == jdebug)
                        {
                            testout[QudPt] = Kernel(chi_tau, chi_t, chi_t - chi_tau);
                        } */
                    }
                    GalerkinMatrix[i + TestDim * j] += LocalMatrix(ii, jj);

                    // Atomic update of the galerkin matrix
                    /* double contribution = LocalMatrix(ii, jj);
                    atomicAdd(&GalerkinMatrix[i + TestDim * j], contribution); */
                }
            }
        }
    }
}