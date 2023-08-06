#include <eigen3/Eigen/Dense>

Eigen::MatrixXd panel_assembly(
    int TrialDim, int TestDim, int NTriangles, int NVertices, int NInteractions,
    int NThreads, const int *I, const int *J, const int *relation,
    const double *W0, const double *X0, int Nq0, const double *W1,
    const double *X1, int Nq1, const double *W2, const double *X2, int Nq2,
    const double *W3, const double *X3, int Nq3, double *shapeDerivative,
    double *GalerkinMatrix, const double *trial_vec, const double *test_vec,
    const int *Elements, const double *Vertices, const double *Normals,
    const double *Areas, const int *Elt2DofTest, const int *Elt2DofTrial,
    int TrialSpace, int TestSpace, int TrialOperator, int TestOperator,
    int NRSFTrial, int NRSFTest) {

  // Looping over all interactions specified by the index pairs I(idx),J(idx)
  for (int InteractionIdx = 0;
       InteractionIdx < NInteractions, ++InteractionIdx) {
    Eigen::Vector3d Ai, Bi, Ci, Aj, Bj, Cj;
    Eigen::MatrixXd Ei(3, 2), Ej(3, 2);

    // Fixed size
    int intersection[3], diffI[3], diffJ[3];

    const double *W = NULL;
    const double *X = NULL;
    int NQudPts = 0;

    // The pair of panels
    int i = I[InteractionIdx], j = J[InteractionIdx];

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

    if (TrialSpace == 0 && TestSpace == 0) {
      // P0 X P0
      for (int ii = 0; ii < NRSFTest; ++ii) {
        double Psix = g_tau;
        for (int jj = 0; jj < NRSFTrial; ++jj) {
          double Psiy = g_t;
          for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {
            Eigen::Vector3d chi_tau =
                Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
            Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[4 * QudPt + 2] +
                                    Ej.col(1) * X[4 * QudPt + 3];

            LocalMatrix(ii, jj) += W[QudPt] * Psix *
                                   SLKernel(chi_tau, chi_t, chi_t - chi_tau) *
                                   Psiy;
          }
          GalerkinMatrix[i + TestDim * j] += LocalMatrix(ii, jj);

          // Atomic update of the galerkin matrix
          // double contribution = LocalMatrix(ii, jj);
          // atomicAdd(&GalerkinMatrix[i + TestDim * j], contribution);
        }
      }
    } else if (TrialSpace == 1 && TestSpace == 1) {
      // P1 X P1

      for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {
        Eigen::Vector3d RSFsX(1 - X[4 * QudPt] - X[4 * QudPt + 1], X[4 * QudPt],
                              X[4 * QudPt + 1]);
        Eigen::Vector3d RSFsY(1 - X[4 * QudPt + 2] - X[4 * QudPt + 3],
                              X[4 * QudPt + 2], X[4 * QudPt + 3]);

        RSFsX *= g_tau;
        RSFsY *= g_t;

        Eigen::Vector3d chi_tau =
            Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
        Eigen::Vector3d chi_t =
            Aj + Ej.col(0) * X[4 * QudPt + 2] + Ej.col(1) * X[4 * QudPt + 3];

        for (int ii = 0; ii < NRSFTest; ++ii) {
          for (int jj = 0; jj < NRSFTrial; ++jj) {
            LocalMatrix(ii, jj) += W[QudPt] * RSFsX(ii) *
                                   SLKernel(chi_tau, chi_t, chi_t - chi_tau) *
                                   RSFsY(jj);
          }
        }
      }
      GalerkinMatrix[EltI[ii] + TestDim * EltJ[jj]] += LocalMatrix(ii, jj);

    } else if (TrialSpace == 1 && TestSpace == 0) {
      // P0 X ntimes(P1) (test X trial)

      for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {

        double Psix = g_tau;
        Eigen::Vector3d RSFsY(1 - X[4 * QudPt + 2] - X[4 * QudPt + 3],
                              X[4 * QudPt + 2], X[4 * QudPt + 3]);

        RSFsY *= g_t;

        Eigen::Vector3d chi_tau =
            Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
        Eigen::Vector3d chi_t =
            Aj + Ej.col(0) * X[4 * QudPt + 2] + Ej.col(1) * X[4 * QudPt + 3];

        for (int ii = 0; ii < NRSFTest; ++ii) {
          for (int jj = 0; jj < NRSFTrial; ++jj) {
            LocalMatrix(ii, jj) +=
                W[QudPt] * Psix *
                SLKernel(chi_tau, chi_t, chi_t - chi_tau).dot(normaly) *
                RSFsY(jj);
          }
        }
      }
      GalerkinMatrix[i + TestDim * EltJ[jj]] += LocalMatrix(ii, jj);
    } else if (TrialSpace == 2 && TestSpace == 2) {
      // RWG X RWG
      for (int ii = 0; ii < NRSFTest; ++ii) {
        int iip1 = (permI[ii] + 1) % 3;
        int iip2 = (iip1 + 1) % 3;

        double fluxI = origEltI[iip1] < origEltI[iip2] ? 1. : -1.;

        for (int jj = 0; jj < NRSFTrial; ++jj) {
          int jjp1 = (permJ[jj] + 1) % 3;
          int jjp2 = (jjp1 + 1) % 3;

          double fluxJ = origEltJ[jjp1] < origEltJ[jjp2] ? 1. : -1.;

          for (int QudPt = 0; QudPt < NQudPts; ++QudPt) {
            /* if (blockIdx.x == 0 && threadIdx.x == 0)
            {
                printf("Qud pt %d\n", QudPt);
            } */

            double RWGX_ref_0 = X[4 * QudPt] - ii % 2;
            double RWGX_ref_1 = X[4 * QudPt + 1] - ii / 2;

            double RWGY_ref_0 = X[4 * QudPt + 2] - jj % 2;
            double RWGY_ref_1 = X[4 * QudPt + 3] - jj / 2;

            Eigen::Vector3d Psix =
                fluxI * (Ei.col(0) * RWGX_ref_0 + Ei.col(1) * RWGX_ref_1);
            Eigen::Vector3d Psiy =
                fluxJ * (Ej.col(0) * RWGY_ref_0 + Ej.col(1) * RWGY_ref_1);

            Eigen::Vector3d chi_tau =
                Ai + Ei.col(0) * X[4 * QudPt] + Ei.col(1) * X[4 * QudPt + 1];
            Eigen::Vector3d chi_t = Aj + Ej.col(0) * X[4 * QudPt + 2] +
                                    Ej.col(1) * X[4 * QudPt + 3];

            LocalMatrix(ii, jj) +=
                W[QudPt] *
                (DLKernel(chi_tau, chi_t, chi_t - chi_tau).cross(Psiy))
                    .dot(Psix);
          }
          GalerkinMatrix[DofsI[permI[ii]] + TestDim * DofsJ[permJ[jj]]] +=
              LocalMatrix(ii, jj);
        }
      }
    }
  }
}