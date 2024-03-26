#include "mex.hpp"
#include "mexAdapter.hpp"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

class MexFunction : public matlab::mex::Function {
public:
  double SLKernel(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX) {
    // return 1;
    return 1. / (4 * M_PI) / YmX.norm();
  }

  Eigen::Vector3d DLKernel(Eigen::Vector3d X, Eigen::Vector3d Y,
                           Eigen::Vector3d YmX) {
    // return 1;
    double znorm = YmX.norm();
    return -YmX / (4 * M_PI) / (znorm * znorm * znorm);
  }

  // Evaluates the block matrix using lowest order BEM spaces (P0,P1)
  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    // Parsing the inputs into matlab data type

    // No. of triangles, vertices and interactions
    matlab::data::TypedArray<int> _NTriangles = std::move(inputs[0]);
    int NTriangles = _NTriangles[0];
    matlab::data::TypedArray<int> _NVertices = std::move(inputs[1]);
    int NVertices = _NVertices[0];
    matlab::data::TypedArray<int> _NInteractions = std::move(inputs[2]);
    int NInteractions = _NInteractions[0];

    // Interaction pair (I,J)
    matlab::data::TypedArray<int> I = std::move(inputs[3]);
    matlab::data::TypedArray<int> J = std::move(inputs[4]);
    matlab::data::TypedArray<int> relation = std::move(inputs[5]);

    // SS Quadrature weights, nodes and number of points
    matlab::data::TypedArray<double> _W0 = std::move(inputs[6]);
    double *W0 = &*_W0.begin();
    matlab::data::TypedArray<double> _X0 = std::move(inputs[7]);
    double *X0 = &*_X0.begin();
    matlab::data::TypedArray<int> _Nq0 = std::move(inputs[8]);
    int Nq0 = _Nq0[0];
    matlab::data::TypedArray<double> _W1 = std::move(inputs[9]);
    double *W1 = &*_W1.begin();
    matlab::data::TypedArray<double> _X1 = std::move(inputs[10]);
    double *X1 = &*_X1.begin();
    matlab::data::TypedArray<int> _Nq1 = std::move(inputs[11]);
    int Nq1 = _Nq1[0];
    matlab::data::TypedArray<double> _W2 = std::move(inputs[12]);
    double *W2 = &*_W2.begin();
    matlab::data::TypedArray<double> _X2 = std::move(inputs[13]);
    double *X2 = &*_X2.begin();
    matlab::data::TypedArray<int> _Nq2 = std::move(inputs[14]);
    int Nq2 = _Nq2[0];
    matlab::data::TypedArray<double> _W3 = std::move(inputs[15]);
    double *W3 = &*_W3.begin();
    matlab::data::TypedArray<double> _X3 = std::move(inputs[16]);
    double *X3 = &*_X3.begin();
    matlab::data::TypedArray<int> _Nq3 = std::move(inputs[17]);
    int Nq3 = _Nq3[0];

    // Elements, vertices, normals and areas
    matlab::data::TypedArray<int> _Elements = std::move(inputs[18]);
    int *Elements = &*_Elements.begin();
    matlab::data::TypedArray<double> _Vertices = std::move(inputs[19]);
    double *Vertices = &*_Vertices.begin();
    matlab::data::TypedArray<double> _Normals = std::move(inputs[20]);
    double *Normals = &*_Normals.begin();
    matlab::data::TypedArray<double> Areas = std::move(inputs[21]);

    // Elt2dof matrices
    matlab::data::TypedArray<int> _Elt2DofTest = std::move(inputs[22]);
    int *Elt2DofTest = &*_Elt2DofTest.begin();
    matlab::data::TypedArray<int> _Elt2DofTrial = std::move(inputs[23]);
    int *Elt2DofTrial = &*_Elt2DofTrial.begin();

    // Permutations of the elements
    matlab::data::TypedArray<int> _PermII = std::move(inputs[24]);
    int *PermII = &*_PermII.begin();
    matlab::data::TypedArray<int> _PermJJ = std::move(inputs[25]);
    int *PermJJ = &*_PermJJ.begin();

    // Galerkin matrices
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(NTriangles, NTriangles);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(NTriangles, NVertices);
    Eigen::MatrixXd W = Eigen::MatrixXd::Zero(NVertices, NVertices);

    // Looping over all interactions specified by the index pairs
    // I(idx),J(idx)
    for (int InteractionIdx = 0; InteractionIdx < NInteractions;
         ++InteractionIdx) {
      // printf("InteractionIdx = %d \n", InteractionIdx);
      Eigen::Vector3d Ai, Bi, Ci, Aj, Bj, Cj;
      Eigen::MatrixXd Ei(3, 2), Ej(3, 2);

      // Fixed size
      int intersection[3], diffI[3], diffJ[3];

      const double *W = NULL;
      const double *X = NULL;
      int NQudPts = 0;

      // The pair of panels
      int i = I[InteractionIdx], j = J[InteractionIdx];

      // std::cout << "i,j = " << i << ", " << j << std::endl;

      double g_tau = 2 * Areas[i], g_t = 2 * Areas[j];

      Eigen::Vector3d normalx(Normals[3 * i], Normals[3 * i + 1],
                              Normals[3 * i + 2]);
      Eigen::Vector3d normaly(Normals[3 * j], Normals[3 * j + 1],
                              Normals[3 * j + 2]);

      int EltI[] = {Elements[3 * i], Elements[3 * i + 1], Elements[3 * i + 2]};

      int EltJ[] = {Elements[3 * j], Elements[3 * j + 1], Elements[3 * j + 2]};

      int origEltI[] = {EltI[0], EltI[1], EltI[2]};
      int origEltJ[] = {EltJ[0], EltJ[1], EltJ[2]};

      int DofsI[] = {Elt2DofTest[3 * i], Elt2DofTest[3 * i + 1],
                     Elt2DofTest[3 * i + 2]};
      int DofsJ[] = {Elt2DofTrial[3 * j], Elt2DofTrial[3 * j + 1],
                     Elt2DofTrial[3 * j + 2]};

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
      for (int k = 0; k < 3; ++k) {
        permI[k] = permII[3 * InteractionIdx + k];
        permJ[k] = permJJ[3 * InteractionIdx + k];
        EltI[k] = origEltI[permI[k]];
        EltJ[k] = origEltJ[permJ[k]];
      }

      // Vertices of element i
      Ai = Eigen::Vector3d(Vertices[3 * EltI[0]], Vertices[3 * EltI[0] + 1],
                           Vertices[3 * EltI[0] + 2]);
      Bi = Eigen::Vector3d(Vertices[3 * EltI[1]], Vertices[3 * EltI[1] + 1],
                           Vertices[3 * EltI[1] + 2]);
      Ci = Eigen::Vector3d(Vertices[3 * EltI[2]], Vertices[3 * EltI[2] + 1],
                           Vertices[3 * EltI[2] + 2]);

      // Vertices of element j
      Aj = Eigen::Vector3d(Vertices[3 * EltJ[0]], Vertices[3 * EltJ[0] + 1],
                           Vertices[3 * EltJ[0] + 2]);
      Bj = Eigen::Vector3d(Vertices[3 * EltJ[1]], Vertices[3 * EltJ[1] + 1],
                           Vertices[3 * EltJ[1] + 2]);
      Cj = Eigen::Vector3d(Vertices[3 * EltJ[2]], Vertices[3 * EltJ[2] + 1],
                           Vertices[3 * EltJ[2] + 2]);

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

      double localvalV = 0;
      Eigen::MatrixXd localMatrixK = Eigen::MatrixXd::Zero(1, 3);
      Eigen::MatrixXd localMatrixW = Eigen::MatrixXd::Zero(3, 3);

      // P0 basis functions
      double Psix_P0 = g_tau;
      double Psiy_P0 = g_t;

      // Computing the local interactions
      for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {
        Eigen::Vector3d chi_tau =
            Ai + Ei * Eigen::Vector2d(X[4 * QudPt], X[4 * QudPt + 1]);
        Eigen::Vector3d chi_t =
            Aj + Ej * Eigen::Vector2d(X[4 * QudPt + 2], X[4 * QudPt + 3]);

        localvalV += Psix_P0 * Psiy_P0 *
                     SLKernel(chi_tau, chi_t, chi_t - chi_tau) * W(QudPt);

        // Trial basis P1
        for (int jj = 0; jj < 3; ++jj) {
          // no. jj rsf
          int Psiy_nxgradP1_0 = (3 - jj) % 3 - 1;
          int Psiy_nxgradP1_1 = jj - 1;

          Eigen::Vector3d Psiy_nxgradP1 =
              g_t * normaly.cross(DCVj * Eigen::Vector2d(Psiy_nxgradP1_0,
                                                         Psiy_nxgradP1_1));

          Eigen::Vector3d RSFsY_P1(1 - X[4 * QudPt + 2] - X[4 * QudPt + 3],
                                   X[4 * QudPt + 2], X[4 * QudPt + 3]);
          RSFsY_P1 *= g_t;

          localMatrixK(0, jj) +=
              W(QudPt) *
              DLKernel(chi_tau, chi_t, chi_t - chi_tau).dot(normaly) *
              RSFsY_P1(jj) * Psix_P0;

          // Test basis P1
          for (int ii = 0; ii < 3; ++ii) {
            // no. ii rsf
            int Psix_nxgradP1_0 = (3 - ii) % 3 - 1;
            int Psix_nxgradP1_1 = ii - 1;

            Eigen::Vector3d Psix_nxgradP1 =
                g_tau * normalx.cross(DCVi * Eigen::Vector2d(Psix_nxgradP1_0,
                                                             Psix_nxgradP1_1));

            Eigen::Vector3d RSFsX_P1(1 - X[4 * QudPt] - X[4 * QudPt + 1],
                                     X[4 * QudPt], [4 * QudPt + 1]);
            RSFsX_P1 *= g_tau;

            localMatrixW(ii, jj) += W(QudPt) *
                                    SLKernel(chi_tau, chi_t, chi_t - chi_tau) *
                                    Psiy_nxgradP1.dot(Psix_nxgradP1);
          }
        }
      }

      // Local to global mapping
      V(i, j) += localvalV;

      for (int jj = 0; jj < 3; ++jj) {
        K(i, EltJ[jj]) += localMatrixK(0, jj);
        for (int ii = 0; ii < 3; ++ii) {
          W(EltI[ii], EltJ[jj]) += localMatrixW(ii, jj);
        }
      }
    }

    // create a MATLAB matrix using factory
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<double> outputMatrix = factory.createArray<double>(
        {NElements + NVertices, NElements + NVertices});

    // Filling up the output Matrix

    // Filling V
    for (int i = 0; i < NElements; ++i) {
      for (int j = 0; j < NElements; ++j) {
        outputMatrix[i][j] = V(i, j);
      }
    }

    // Filling K and K'
    for (int i = 0; i < NElements; ++i) {
      for (int j = 0; j < NVertices; ++j) {
        outputMatrix[i][NElements + j] = K(i, j);
        outputMatrix[NElements + i][j] = K(j, i);
      }
    }

    // Filling W
    for (int i = 0; i < NVertices; ++i) {
      for (int j = 0; j < NVertices; ++j) {
        outputMatrix[NElements + i][NElements + j] = W(i, j);
      }
    }

    outputs[0] = outputMatrix;
  }
};
