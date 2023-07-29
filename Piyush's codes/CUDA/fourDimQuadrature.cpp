#include <cassert>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <stdlib.h>

using namespace std;

Eigen::VectorXd phi(int maxDegree, Eigen::VectorXd &weights,
                    Eigen::VectorXd &X1, Eigen::VectorXd &X2,
                    Eigen::VectorXd &X3, Eigen::VectorXd &X4) {
  int Nqud = weights.size();
  assert(X1.size() == Nqud);
  assert(X2.size() == Nqud);
  assert(X3.size() == Nqud);
  assert(X4.size() == Nqud);
  assert(maxDegree * maxDegree * maxDegree * maxDegree == 5 * Nqud);

  Eigen::VectorXd out =
      Eigen::VectorXd::Zero(maxDegree * maxDegree * maxDegree * maxDegree);

  // Evaluate all the polynomials
  for (int i = 0; i < maxDegree; ++i) {       // exponent of x1
    for (int j = 0; j < maxDegree; ++j) {     // exponent of x2
      for (int k = 0; k < maxDegree; ++k) {   // exponent of x3
        for (int l = 0; l < maxDegree; ++l) { // exponent of x4
          // Evaluating quadrature
          for (int idx = 0; idx < Nqud; ++idx) {
            out(maxDegree * maxDegree * maxDegree * i +
                maxDegree * maxDegree * j + maxDegree * k + l) +=
                weights[idx] * pow(X1(idx), i) * pow(X2(idx), j) *
                pow(X3(idx), k) * pow(X4(idx), l);
          }
          // subtracting the RHS
          out(maxDegree * maxDegree * maxDegree * i +
              maxDegree * maxDegree * j + maxDegree * k + l) -=
              1. / (i + 1) * 1. / (j + 1) * 1. / (k + 1) * 1. / (l + 1);
        }
      }
    }
  }
  return out;
}

Eigen::MatrixXd Dphi(int maxDegree, Eigen::VectorXd &weights,
                     Eigen::VectorXd &X1, Eigen::VectorXd &X2,
                     Eigen::VectorXd &X3, Eigen::VectorXd &X4) {
  int Nqud = weights.size();
  Eigen::MatrixXd out =
      Eigen::MatrixXd::Zero(maxDegree * maxDegree * maxDegree * maxDegree,
                            maxDegree * maxDegree * maxDegree * maxDegree);

  // The first four loops give the row idx
  for (int i = 0; i < maxDegree; ++i) {       // exponent of x1
    for (int j = 0; j < maxDegree; ++j) {     // exponent of x2
      for (int k = 0; k < maxDegree; ++k) {   // exponent of x3
        for (int l = 0; l < maxDegree; ++l) { // exponent of x4
          // Evaluating quadrature
          for (int idx = 0; idx < Nqud; ++idx) {
            // Derivative with respect to weights
            out(maxDegree * maxDegree * maxDegree * i +
                    maxDegree * maxDegree * j + maxDegree * k + l,
                idx) = pow(X1(idx), i) * pow(X2(idx), j) * pow(X3(idx), k) *
                       pow(X4(idx), l);

            // Derivative with respect to X1
            out(maxDegree * maxDegree * maxDegree * i +
                    maxDegree * maxDegree * j + maxDegree * k + l,
                idx + Nqud) = i * weights[idx] * pow(X1(idx), i - 1) *
                              pow(X2(idx), j) * pow(X3(idx), k) *
                              pow(X4(idx), l);

            // Derivative with respect to X2
            out(maxDegree * maxDegree * maxDegree * i +
                    maxDegree * maxDegree * j + maxDegree * k + l,
                idx + 2 * Nqud) = j * weights[idx] * pow(X1(idx), i) *
                                  pow(X2(idx), j - 1) * pow(X3(idx), k) *
                                  pow(X4(idx), l);

            // Derivative with respect to X3
            out(maxDegree * maxDegree * maxDegree * i +
                    maxDegree * maxDegree * j + maxDegree * k + l,
                idx + 3 * Nqud) = k * weights[idx] * pow(X1(idx), i) *
                                  pow(X2(idx), j) * pow(X3(idx), k - 1) *
                                  pow(X4(idx), l);

            // Derivative with respect to X4
            out(maxDegree * maxDegree * maxDegree * i +
                    maxDegree * maxDegree * j + maxDegree * k + l,
                idx + 4 * Nqud) = l * weights[idx] * pow(X1(idx), i) *
                                  pow(X2(idx), j) * pow(X3(idx), k) *
                                  pow(X4(idx), l - 1);
          }
        }
      }
    }
  }
  return out;
}

void newtonIteration(Eigen::VectorXd &x0, int maxDegree, int maxit = 100,
                     double atol = 1e-10, double rtol = 1e-12) {

  int N = x0.size();
  int Nbyfive = N / 5;
  cout << "Nbyfive = " << Nbyfive << endl;

  auto phiHelper = [&](Eigen::VectorXd x) {
    Eigen::VectorXd wts = x.segment(0, Nbyfive);
    Eigen::VectorXd x1 = x.segment(Nbyfive, Nbyfive);
    Eigen::VectorXd x2 = x.segment(2 * Nbyfive, Nbyfive);
    Eigen::VectorXd x3 = x.segment(3 * Nbyfive, Nbyfive);
    Eigen::VectorXd x4 = x.segment(4 * Nbyfive, Nbyfive);
    return phi(maxDegree, wts, x1, x2, x3, x4);
  };

  auto DphiHelper = [&](Eigen::VectorXd x) {
    Eigen::VectorXd wts = x.segment(0, Nbyfive);
    Eigen::VectorXd x1 = x.segment(Nbyfive, Nbyfive);
    Eigen::VectorXd x2 = x.segment(2 * Nbyfive, Nbyfive);
    Eigen::VectorXd x3 = x.segment(3 * Nbyfive, Nbyfive);
    Eigen::VectorXd x4 = x.segment(4 * Nbyfive, Nbyfive);
    return Dphi(maxDegree, wts, x1, x2, x3, x4);
  };

  Eigen::MatrixXd Jacobian = DphiHelper(x0);
  // cout << "Jacobian \n" << Jacobian << endl;
  Eigen::VectorXd Vec = phiHelper(x0);
  Eigen::VectorXd correction = Jacobian.lu().solve(Vec);
  cout << "Jacobian determinant = " << Jacobian.determinant()
       << " , vec norm = " << Vec.norm() << endl;
  int it = 0;
  cout << "correction.norm() = " << correction.norm() << endl;
  while (correction.norm() > atol || it < maxit) {
    cout << "iteration no. " << it << " correction norm = " << correction.norm()
         << endl;
    x0 -= correction;
    Eigen::VectorXd correction = DphiHelper(x0).lu().solve(phiHelper(x0));
    it++;
  }
}

int main() {
  int Nsystem = 625;
  int maxDegree = 5;
  srand(99);
  Eigen::VectorXd wtsAndPts = Eigen::VectorXd::Random(Nsystem);
  newtonIteration(wtsAndPts, maxDegree);
  return 0;
}