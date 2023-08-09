#include "mex.hpp"
#include "mexAdapter.hpp"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <vector>

class MexFunction : public matlab::mex::Function {
public:
  Eigen::Vector3d Vel(int a, int b, int c, int alpha, Eigen::Vector3d X) {
    Eigen::Vector3d out = Eigen::Vector3d::Zero(3);
    out(alpha) = cos(a * X(0)) * cos(b * X(1)) * cos(c * X(2));
    return out;
  }

  Eigen::Matrix3d DVel(int a, int b, int c, int alpha, Eigen::Vector3d X) {
    Eigen::Matrix3d out = Eigen::Matrix3d::Zero(3, 3);
    out.row(alpha) =
        -Eigen::Vector3d(a * sin(a * X(0)) * cos(b * X(1)) * cos(c * X(2)),
                         b * cos(a * X(0)) * sin(b * X(1)) * cos(c * X(2)),
                         c * cos(a * X(0)) * cos(b * X(1)) * sin(c * X(2)));
    return out;
  }

  double KernelA1(int a, int b, int c, int alpha, Eigen::Vector3d X,
                  Eigen::Vector3d Y, Eigen::Vector3d YmX) {
    double znorm = YmX.norm();
    return YmX.dot(Vel(a, b, c, alpha, X) - Vel(a, b, c, alpha, Y)) /
           (4 * M_PI) / (znorm * znorm * znorm);
  }

  double KernelA2(Eigen::Vector3d X, Eigen::Vector3d Y, Eigen::Vector3d YmX) {
    return 1. / (4 * M_PI) / YmX.norm();
  }

  Eigen::Vector3d KernelC1(Eigen::Vector3d X, Eigen::Vector3d Y,
                           Eigen::Vector3d YmX) {
    // return 1;
    double znorm = YmX.norm();
    return YmX / (4 * M_PI) / (znorm * znorm * znorm);
  }

  Eigen::Vector3d KernelC3(int a, int b, int c, int alpha, Eigen::Vector3d X,
                           Eigen::Vector3d Y, Eigen::Vector3d YmX) {
    // return 1;
    double znorm = YmX.norm();
    return -3. / (4 * M_PI) * YmX *
               YmX.dot(Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) /
               (znorm * znorm * znorm * znorm * znorm) //
           +
           1. / (4 * M_PI) * (Vel(a, b, c, alpha, Y) - Vel(a, b, c, alpha, X)) /
               (znorm * znorm * znorm);
  }

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

  void IntersectionDiff(int *EltI, int *EltJ, int intersection[], int diffI[],
                        int diffJ[]) {
    bool EltITracker[] = {false, false, false};
    bool EltJTracker[] = {false, false, false};

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (EltI[i] == EltJ[j]) {
          EltITracker[i] = true;
          EltJTracker[j] = true;
        }
      }
    }

    int interidx = 0, diffiidx = 0, diffjidx = 0;

    for (int i = 0; i < 3; ++i) {
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

  void operator()(matlab::mex::ArgumentList outputs,
                  matlab::mex::ArgumentList inputs) {
    // Parsing the inputs into matlab data type
    matlab::data::TypedArray<int> _TrialDim = std::move(inputs[0]);
    int TrialDim = _TrialDim[0];
    matlab::data::TypedArray<int> _TestDim = std::move(inputs[1]);
    int TestDim = _TestDim[0];
    matlab::data::TypedArray<int> _NTriangles = std::move(inputs[2]);
    int NTriangles = _NTriangles[0];
    matlab::data::TypedArray<int> _NVertices = std::move(inputs[3]);
    int NVertices = _NVertices[0];
    matlab::data::TypedArray<int> _NInteractions = std::move(inputs[4]);
    int NInteractions = _NInteractions[0];
    matlab::data::TypedArray<int> I = std::move(inputs[5]);
    matlab::data::TypedArray<int> J = std::move(inputs[6]);
    matlab::data::TypedArray<int> relation = std::move(inputs[7]);
    matlab::data::TypedArray<double> _W0 = std::move(inputs[8]);
    double *W0 = &*_W0.begin();
    matlab::data::TypedArray<double> _X0 = std::move(inputs[9]);
    double *X0 = &*_X0.begin();
    matlab::data::TypedArray<int> _Nq0 = std::move(inputs[10]);
    int Nq0 = _Nq0[0];
    matlab::data::TypedArray<double> _W1 = std::move(inputs[11]);
    double *W1 = &*_W1.begin();
    matlab::data::TypedArray<double> _X1 = std::move(inputs[12]);
    double *X1 = &*_X1.begin();
    matlab::data::TypedArray<int> _Nq1 = std::move(inputs[13]);
    int Nq1 = _Nq1[0];
    matlab::data::TypedArray<double> _W2 = std::move(inputs[14]);
    double *W2 = &*_W2.begin();
    matlab::data::TypedArray<double> _X2 = std::move(inputs[15]);
    double *X2 = &*_X2.begin();
    matlab::data::TypedArray<int> _Nq2 = std::move(inputs[16]);
    int Nq2 = _Nq2[0];
    matlab::data::TypedArray<double> _W3 = std::move(inputs[17]);
    double *W3 = &*_W3.begin();
    matlab::data::TypedArray<double> _X3 = std::move(inputs[18]);
    double *X3 = &*_X3.begin();
    matlab::data::TypedArray<int> _Nq3 = std::move(inputs[19]);
    int Nq3 = _Nq3[0];
    matlab::data::TypedArray<double> trial_vec = std::move(inputs[20]);
    matlab::data::TypedArray<double> test_vec = std::move(inputs[21]);
    matlab::data::TypedArray<int> _Elements = std::move(inputs[22]);
    int *Elements = &*_Elements.begin();
    matlab::data::TypedArray<double> _Vertices = std::move(inputs[23]);
    double *Vertices = &*_Vertices.begin();
    matlab::data::TypedArray<double> _Normals = std::move(inputs[24]);
    double *Normals = &*_Normals.begin();
    matlab::data::TypedArray<double> Areas = std::move(inputs[25]);
    matlab::data::TypedArray<int> _Elt2DofTest = std::move(inputs[26]);
    int *Elt2DofTest = &*_Elt2DofTest.begin();
    matlab::data::TypedArray<int> _Elt2DofTrial = std::move(inputs[27]);
    int *Elt2DofTrial = &*_Elt2DofTrial.begin();
    matlab::data::TypedArray<int> _TrialSpace = std::move(inputs[28]);
    int TrialSpace = _TrialSpace[0];
    matlab::data::TypedArray<int> _TestSpace = std::move(inputs[29]);
    int TestSpace = _TestSpace[0];
    matlab::data::TypedArray<int> _TrialOperator = std::move(inputs[30]);
    int TrialOperator = _TrialOperator[0];
    matlab::data::TypedArray<int> _TestOperator = std::move(inputs[31]);
    int TestOperator = _TestOperator[0];
    matlab::data::TypedArray<int> _NRSFTrial = std::move(inputs[32]);
    int NRSFTrial = _NRSFTrial[0];
    matlab::data::TypedArray<int> _NRSFTest = std::move(inputs[33]);
    int NRSFTest = _NRSFTest[0];
    matlab::data::TypedArray<double> TdA = std::move(inputs[34]);
    matlab::data::TypedArray<double> TnA = std::move(inputs[35]);
    matlab::data::TypedArray<double> Mxn_coeffs = std::move(inputs[36]);
    matlab::data::TypedArray<int> _abc_alpha = std::move(inputs[37]);
    int *abc_alpha = &*_abc_alpha.begin();
    matlab::data::TypedArray<int> _Nabc_alpha = std::move(inputs[38]);
    int Nabc_alpha = _Nabc_alpha[0];
    matlab::data::TypedArray<double> _mu = std::move(inputs[39]);
    double mu = _mu[0];
    matlab::data::TypedArray<double> _mu0 = std::move(inputs[40]);
    double mu0 = _mu0[0];

    // Vector storing the output shape derivatives
    std::vector<double> localShapeDerivatives(Nabc_alpha, 0);
    // Looping over all interactions specified by the index pairs I(idx),J(idx)
    for (int InteractionIdx = 0; InteractionIdx < NInteractions;
         ++InteractionIdx) {
      printf("Interaction No. %d \n", InteractionIdx);
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

      if (relation[InteractionIdx] == 0) // No interaction
      {
        // Computing Quadrature
        W = W0;
        X = X0;
        NQudPts = Nq0;
      } else if (relation[InteractionIdx] == 1) // Common vertex
      {
        IntersectionDiff(EltI, EltJ, intersection, diffI, diffJ);

        for (int l = 0; l < 3; ++l) {
          // Permutation for I
          if (EltI[l] == intersection[0]) {
            permI[0] = l;
          } else if (EltI[l] == diffI[0]) {
            permI[1] = l;
          } else if (EltI[l] == diffI[1]) {
            permI[2] = l;
          }

          // Permutation for J
          if (EltJ[l] == intersection[0]) {
            permJ[0] = l;
          } else if (EltJ[l] == diffJ[0]) {
            permJ[1] = l;
          } else if (EltJ[l] == diffJ[1]) {
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
      } else if (relation[InteractionIdx] == 2) // Common edge
      {
        IntersectionDiff(EltI, EltJ, intersection, diffI, diffJ);

        for (int l = 0; l < 3; ++l) {
          // Permutation for I
          if (EltI[l] == intersection[0]) {
            permI[0] = l;
          } else if (EltI[l] == intersection[1]) {
            permI[1] = l;
          } else if (EltI[l] == diffI[0]) {
            permI[2] = l;
          }

          // Permutation for J
          if (EltJ[l] == intersection[0]) {
            permJ[0] = l;
          } else if (EltJ[l] == intersection[1]) {
            permJ[1] = l;
          } else if (EltJ[l] == diffJ[0]) {
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
      } else // Identical panels, case 3
      {
        // Computing Quadrature
        W = W3;
        X = X3;
        NQudPts = Nq3;
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

      Eigen::MatrixXd LocalMatrix = Eigen::MatrixXd::Zero(NRSFTest, NRSFTrial);

      for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {
        Eigen::Vector3d chi_tau =
            Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
        Eigen::Vector3d chi_t =
            Aj + Ej.col(0) * X[4 * QudPt + 2] + Ej.col(1) * X[4 * QudPt + 3];

        // Looping over different fields
        for (int fieldIdx = 0; fieldIdx < Nabc_alpha; ++fieldIdx) {
          double localval = 0;
          for (int ii = 0; ii < 3; ++ii) {
            int iip1 = (permI[ii] + 1) % 3;
            int iip2 = (iip1 + 1) % 3;

            double fluxI = origEltI[iip1] < origEltI[iip2] ? 1. : -1.;
            for (int jj = 0; jj < 3; ++jj) {
              int jjp1 = (permJ[jj] + 1) % 3;
              int jjp2 = (jjp1 + 1) % 3;

              double fluxJ = origEltJ[jjp1] < origEltJ[jjp2] ? 1. : -1.;

              double RWGX_ref_0 = X[4 * QudPt] - ii % 2;
              double RWGX_ref_1 = X[4 * QudPt + 1] - ii / 2;

              double RWGY_ref_0 = X[4 * QudPt + 2] - jj % 2;
              double RWGY_ref_1 = X[4 * QudPt + 3] - jj / 2;

              // RWG elements
              Eigen::Vector3d Psix =
                  fluxI * (Ei.col(0) * RWGX_ref_0 + Ei.col(1) * RWGX_ref_1);
              Eigen::Vector3d Psiy =
                  fluxJ * (Ej.col(0) * RWGY_ref_0 + Ej.col(1) * RWGY_ref_1);

              // A1 with RWG . RWG
              double LocalMatrixA1_ii_jj =
                  W[QudPt] //
                  *
                  KernelA1(abc_alpha[4 * fieldIdx], abc_alpha[4 * fieldIdx + 1],
                           abc_alpha[4 * fieldIdx + 2],
                           abc_alpha[4 * fieldIdx + 3], chi_tau, chi_t,
                           chi_t - chi_tau) *
                  Psiy.dot(Psix);

              // A2 with DVelRWG(y) . RWG(x)
              double LocalMatrixA2_ii_jj =
                  W[QudPt] //
                  * KernelA2(chi_tau, chi_t, chi_t - chi_tau) *
                  (DVel(abc_alpha[4 * fieldIdx],     //
                        abc_alpha[4 * fieldIdx + 1], //
                        abc_alpha[4 * fieldIdx + 2], //
                        abc_alpha[4 * fieldIdx + 3], //
                        chi_t) *
                   Psiy)
                      .dot(Psix);

              // C1 with DVelRWG(y) X RWG(X)
              double LocalMatrixC1_ii_jj =
                  W[QudPt] * //
                  (KernelC1(chi_tau, chi_t, chi_t - chi_tau)
                       .cross(DVel(abc_alpha[4 * fieldIdx],     //
                                   abc_alpha[4 * fieldIdx + 1], //
                                   abc_alpha[4 * fieldIdx + 2], //
                                   abc_alpha[4 * fieldIdx + 3], chi_t) *
                              Psiy))
                      .dot(Psix);

              // C3 with RWG X RWG
              double LocalMatrixC3_ii_jj =
                  W[QudPt]                                 //
                  * (KernelC3(abc_alpha[4 * fieldIdx],     //
                              abc_alpha[4 * fieldIdx + 1], //
                              abc_alpha[4 * fieldIdx + 2], //
                              abc_alpha[4 * fieldIdx + 3], //
                              chi_tau, chi_t, chi_t - chi_tau)
                         .cross(Psiy))
                        .dot(Psix);

              // N with RWG.div X RWG.div, kernelN = kernelA1
              double LocalMatrixN_ii_jj =
                  W[QudPt]                                //
                  * KernelA1(abc_alpha[4 * fieldIdx],     //
                             abc_alpha[4 * fieldIdx + 1], //
                             abc_alpha[4 * fieldIdx + 2], //
                             abc_alpha[4 * fieldIdx + 3], //
                             chi_tau, chi_t, chi_t - chi_tau) *
                  4 * fluxI * fluxJ;

              // Combining all the atomic adds above into a single one
              localShapeDerivatives[fieldIdx] +=
                  1 / (2 * mu0) *
                      ((1 + mu / mu0) *
                           (TnA[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj *
                                TnA[DofsJ[permJ[jj]]] // A1
                            + 2 * TnA[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj *
                                  TnA[DofsJ[permJ[jj]]]) // A2
                       - 4 * (TnA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj *
                                  TdA[DofsJ[permJ[jj]]] // C1
                              + TdA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj *
                                    TnA[DofsJ[permJ[jj]]] // C2
                              + TnA[DofsI[permI[ii]]] * LocalMatrixC3_ii_jj *
                                    TdA[DofsJ[permJ[jj]]]) // C3
                       + (1 + mu0 / mu) * -TdA[DofsI[permI[ii]]] *
                             LocalMatrixN_ii_jj * TdA[DofsJ[permJ[jj]]]) // N
                  +
                  mu / 2 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj *
                      Mxn_coeffs[DofsJ[permJ[jj]]] // red_remaining
                  + mu * TnA[DofsI[permI[ii]]] * LocalMatrixA1_ii_jj *
                        Mxn_coeffs[DofsJ[permJ[jj]]] // red_l_M
                  + mu * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj *
                        Mxn_coeffs[DofsJ[permJ[jj]]] // blue_remaining
                  + mu * (Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj *
                              TnA[DofsJ[permJ[jj]]] + // blue_l_M
                          TnA[DofsI[permI[ii]]] * LocalMatrixA2_ii_jj *
                              Mxn_coeffs[DofsJ[permJ[jj]]]) +
                  -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj *
                      TdA[DofsJ[permJ[jj]]] // MC1
                  + -mu0 * TdA[DofsI[permI[ii]]] * LocalMatrixC1_ii_jj *
                        Mxn_coeffs[DofsJ[permJ[jj]]] // MC2
                  + -mu0 * Mxn_coeffs[DofsI[permI[ii]]] * LocalMatrixC3_ii_jj *
                        TdA[DofsJ[permJ[jj]]]; // MC3
            }
          }
          // Combining all the atomic adds above into a single one
          // atomicAdd(&shapeDerivative[fieldIdx], localval);
        }
      }
    }

    // create a MATLAB matrix using factory
    matlab::data::ArrayFactory factory;
    matlab::data::TypedArray<double> output =
        factory.createArray<double>({Nabc_alpha});

    // Filling up the output Matrix
    for (int i = 0; i < Nabc_alpha; ++i) {
      output[i] = localShapeDerivatives[i];
    }

    outputs[0] = output;
  }
};
