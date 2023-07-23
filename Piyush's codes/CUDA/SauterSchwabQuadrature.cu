#include <cuda_runtime.h>
#include <stdio.h>
#include <cuda.h>
#include <math.h>
#include <thrust/device_vector.h>
#include <thrust/set_operations.h>
#include <thrust/sort.h>

#include <cassert>
#include <iostream>
#include <memory>
#include <eigen3/Eigen/Dense>

__device__ double Kernel(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX)
{
    return 1. / (4 * M_PI) / YmX.norm();
}

// Returns Intersection, DiffI, DiffJ
__device__ thrust::tuple<thrust::device_vector<int>, thrust::device_vector<int>, thrust::device_vector<int>> ComputeIntersectionDiff(int *EltI, int *EltJ)
{
    thrust::device_vector<bool> EltITracker(3, false), EltJTracker(3, false);
    thrust::device_vector<int> intersection{}, diffI{}, diffJ{};
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

    for (int i = 0; i < 3; ++i)
    {
        if (EltITracker[i])
            intersection.push_back(EltI[i]);
        else
            diffI.push_back(EltI[i]);

        if (!EltJTracker[i])
            diffJ.push_back(EltJ[i]);
    }

    // Determine diffs
    thrust::tuple<thrust::device_vector<int>, thrust::device_vector<int>, thrust::device_vector<int>> output(intersection, diffI, diffJ);
    return output;
}

// Defining the kernel for SS computation
/*
 * inputs
 *
 * NTriangles : Number of elements
 */
__global__ void computeShapeDerivative(int TrialDim, int TestDim, int NTriangles,
                                       int NThreads, int *I, int *J, int *relation,
                                       double *W0, double **X0, int Nq0,
                                       double *W1, double **X1, int Nq1,
                                       double *W2, double **X2, int Nq2,
                                       double *W3, double **X3, int Nq3,
                                       double *shapeDerivative,
                                       double **GalerkinMatrix,
                                       double *trial_vec, double *test_vec,
                                       int **Elements, double **Vertices, double **Normals, double *Areas,
                                       int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
                                       int NRSFTrial, int NRSFTest)
{
    int ThreadID = blockIdx.x * blockDim.x + threadIdx.x;

    // Size of the matrix
    // int NInteractions = TrialDim * TestDim;

    // Number of element interaction
    int NInteractions = NTriangles * NTriangles;

    // Matrix size is NTriangles^2. Each entry is assigned to a thread
    // InteractionsPerThread gives the max no.of element interactions assigned
    // to a thread
    int InteractionsPerThread = ceil(double(NInteractions) / double(NThreads));

    // Stores the contribution to the shape derivative from the
    // local matrix that will be computed in this thread
    double localval = 0;

    // Looping over all assigned interactions
    for (int idx = 0; idx < InteractionsPerThread; ++idx)
    {
        // The interaction number
        // int InteractionIdx = ThreadID + NThreads * idx;
        int InteractionIdx = ThreadID * InteractionsPerThread + idx;

        // Preparing variables
        Eigen::Vector3d Ai, Bi, Ci, ABCi, permi, Aj, Bj, Cj, ABCj, permj;
        Eigen::MatrixXd Ei(3, 2), Ej(3, 2);

        double *W = NULL;
        double **X = NULL;
        int NQudPts = 0;

        // Make sure that the last thread stays in limit
        if (InteractionIdx >= NInteractions)
            break;

        // The pair of panels
        int i = I[InteractionIdx], j = J[InteractionIdx];

        double g_tau = 2 * Areas[i], g_t = 2 * Areas[j];

        if (relation[InteractionIdx] == 0) // No interaction
        {
            // Vertices of element i
            Ai = Eigen::Vector3d(Vertices[Elements[i][0]]);
            Bi = Eigen::Vector3d(Vertices[Elements[i][1]]);
            Ci = Eigen::Vector3d(Vertices[Elements[i][2]]);

            // Vertices of element j
            Aj = Eigen::Vector3d(Vertices[Elements[j][0]]);
            Bj = Eigen::Vector3d(Vertices[Elements[j][1]]);
            Cj = Eigen::Vector3d(Vertices[Elements[j][2]]);

            // Computing Quadrature
            W = W0;
            X = X0;
            NQudPts = Nq0;
        }
        else if (relation[InteractionIdx] == 1) // Common vertex
        {
            thrust::device_vector<int> intersection;
            thrust::device_vector<int> diffI;
            thrust::device_vector<int> diffJ;
            thrust::tie(intersection, diffI, diffJ) = ComputeIntersectionDiff(Elements[i], Elements[j]);

            // Vertices of element i
            Ai = Eigen::Vector3d(Vertices[Elements[intersection[0]][0]]);
            Bi = Eigen::Vector3d(Vertices[Elements[diffI[0]][1]]);
            Ci = Eigen::Vector3d(Vertices[Elements[diffI[1]][2]]);

            // Vertices of element j
            Aj = Eigen::Vector3d(Vertices[Elements[intersection[0]][0]]);
            Bj = Eigen::Vector3d(Vertices[Elements[diffJ[0]][1]]);
            Cj = Eigen::Vector3d(Vertices[Elements[diffJ[1]][2]]);

            // Computing Quadrature
            W = W1;
            X = X1;
            NQudPts = Nq1;
        }
        else if (relation[InteractionIdx] == 2) // Common edge
        {
            thrust::device_vector<int> intersection;
            thrust::device_vector<int> diffI;
            thrust::device_vector<int> diffJ;
            thrust::tie(intersection, diffI, diffJ) = ComputeIntersectionDiff(Elements[i], Elements[j]);

            // Vertices of element i
            Ai = Eigen::Vector3d(Vertices[Elements[intersection[0]][0]]);
            Bi = Eigen::Vector3d(Vertices[Elements[intersection[1]][1]]);
            Ci = Eigen::Vector3d(Vertices[Elements[diffI[0]][2]]);

            // Vertices of element j
            Aj = Eigen::Vector3d(Vertices[Elements[intersection[0]][0]]);
            Bj = Eigen::Vector3d(Vertices[Elements[intersection[1]][1]]);
            Cj = Eigen::Vector3d(Vertices[Elements[diffJ[0]][2]]);

            // Computing Quadrature
            W = W2;
            X = X2;
            NQudPts = Nq2;
        }
        else // Identical panels, case 3
        {
            // Vertices of element i
            Ai = Eigen::Vector3d(Vertices[Elements[i][0]]);
            Bi = Eigen::Vector3d(Vertices[Elements[i][1]]);
            Ci = Eigen::Vector3d(Vertices[Elements[i][2]]);

            // Vertices of element j
            Aj = Eigen::Vector3d(Vertices[Elements[j][0]]);
            Bj = Eigen::Vector3d(Vertices[Elements[j][1]]);
            Cj = Eigen::Vector3d(Vertices[Elements[j][2]]);

            // Computing Quadrature
            W = W3;
            X = X3;
            NQudPts = Nq3;
        }

        // Jacobian Matrices

        Ei.col(0) = Bi - Ai;
        Ei.col(1) = Ci - Ai;

        Ej.col(0) = Bj - Aj;
        Ej.col(1) = Cj - Aj;

        Eigen::MatrixXd Dxyi = (Ei.transpose() * Ei).inverse();
        Eigen::MatrixXd Dxyj = (Ej.transpose() * Ej).inverse();
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
                        Eigen::Vector3d chi_tau = Ai + Ei.col(0) * X[QudPt][0] + Ei.col(1) * X[QudPt][1];
                        Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[QudPt][3] + Ej.col(1) * X[QudPt][4];
                        LocalMatrix(ii, jj) += W[QudPt] * Psix * Kernel(chi_tau, chi_t, chi_t - chi_tau) * Psiy;
                    }
                    GalerkinMatrix[i][j] += LocalMatrix(ii, jj);
                }
            }
        }
    }
}

// Link with matlab

int main()
{
    // Indices received from matlab are ordered based on interactions

    return 0;
}