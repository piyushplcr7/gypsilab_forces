#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

Eigen::Vector2d kite(double t) {
  return Eigen::Vector2d(0.3+.35 * std::cos(t) + .1625 * std::cos(2 * t),
                         0.5+.35 * std::sin(t));
}

Eigen::Vector2d dkite(double t) {
  return Eigen::Vector2d(-.35 * std::sin(t) - 2 * .1625 * std::sin(2 * t),
                         .35 * std::cos(t));
}

Eigen::VectorXd get_kite_params(unsigned N) {
  // Calculating the length of the kite
  unsigned N_length = 500; // No. of points used in the calculation
  Eigen::VectorXd pts_length = Eigen::VectorXd::LinSpaced(N_length,0,2*M_PI);
  double L = 0;
  for (unsigned i = 0 ; i < N_length-1 ; ++i)
    L += (kite(pts_length(i)) - kite(pts_length(i+1))).norm();

  std::cout << "found length of the kite: " << L << std::endl;
  // Solving the equation for Phi using explicit timestepping
  unsigned k = 20; // multiplicity factor?
  double h = L/N/k; // step size
  Eigen::VectorXd phi_full = Eigen::VectorXd::Constant(N*k,0);
  Eigen::VectorXd phi = Eigen::VectorXd::Constant(N,0);
  for (unsigned i = 1 ; i < N*k ; ++i)
    phi_full(i) = phi_full(i-1) + h /( dkite(phi_full(i-1)) ).norm();

  for (unsigned i = 0 ; i < N ; ++i)
    phi(i) = phi_full(i*k);

  return phi;
}

std::string add_point(Eigen::Vector2d pt, double h, unsigned idx) {
  std::string lhs = "Point(" + std::to_string(idx) + ") = ";
  std::string rhs = "{" + std::to_string(pt(0)) + "," + std::to_string(pt(1)) +
                    ",0," + std::to_string(h) + "};";
  return lhs + rhs;
}

void add_lines(std::ofstream &out, unsigned N_lines, unsigned pt_start_idx,
               unsigned line_start_idx) {
  std::string new_entity_command = "//+";
  for (unsigned i = 0; i < N_lines; ++i) {
    out << new_entity_command << std::endl;
    std::string lhs = "Line(" + std::to_string(i + line_start_idx) + ") = ";
    std::string rhs = "{" + std::to_string(pt_start_idx + i) + "," +
                      std::to_string(pt_start_idx + (i + 1) % N_lines) + "};";
    out << lhs + rhs << std::endl;
  }
}

void add_curve_loop(std::ofstream& out, unsigned idx, unsigned start, unsigned N) {
  out << "//+" << std::endl;
  std::string lhs = "Curve Loop(" + std::to_string(idx) + ") = ";
  std::string rhs = "{";
  for (unsigned i = 0 ; i < N ; ++i) {
    rhs += std::to_string(i+start);
    if (i != N-1)
      rhs += ",";
  }
  rhs += "};";
  out << lhs + rhs << std::endl;
}

void add_physical_curve(std::ofstream& out, unsigned idx, unsigned start, unsigned N, std::string tag) {
  out << "//+" << std::endl;
  std::string lhs = "Physical Curve(\"" + tag + "\"," + std::to_string(idx) + ") = ";
  std::string rhs = "{";
  for (unsigned i = 0 ; i < N ; ++i) {
    rhs += std::to_string(i+start);
    if (i != N-1)
      rhs += ",";
  }
  rhs += "};";
  out << lhs + rhs << std::endl;
}

int main() {
  // Creating the geometry file
  std::string filename = "kite_sq.geo";
  std::ofstream out(filename);
  std::string new_entity_command = "//+";

  // Mesh size parameters for gmsh
  double h_sq = 1;
  double h_kite = 1;

  // Edge length of the square
  double s = 4;
  // No. of points on one edge of the square (excluding the corner points)
  unsigned M = MM;
  // Working size
  unsigned N = M+2;
  Eigen::VectorXd var_coords = Eigen::VectorXd::LinSpaced(N,-s/2,s/2);
  Eigen::VectorXd const_coords = Eigen::VectorXd::Constant(N,s/2);

  // Right edge of the square
  Eigen::MatrixXd pts_r(N,2);
  pts_r << const_coords, var_coords;

  // Adding the points, starting from number 1 here and adding N points
  for (unsigned i = 0 ; i < N ; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(pts_r.row(i),h_sq,i+1) << std::endl;
  }

  // Top edge of the square
  Eigen::MatrixXd pts_t(N,2);
  pts_t << var_coords.reverse(),const_coords;

  // Adding the points, starting from number N+1 here and adding N-1 points
  for (unsigned i = 1 ; i < N ; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(pts_t.row(i),h_sq,i+N) << std::endl;
  }

  // Left edge of the square
  Eigen::MatrixXd pts_l(N,2);
  pts_l << -const_coords,var_coords.reverse();

  // Adding the points, starting from number 2N here and adding N-1 points
  for (unsigned i = 1 ; i < N ; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(pts_l.row(i),h_sq,i+2*N-1) << std::endl;
  }

  // Bottom edge of the square
  Eigen::MatrixXd pts_b(N,2);
  pts_b << var_coords,-const_coords;

  // Adding the points, starting from number 3N-1 here and adding N-2 points
  for (unsigned i = 1 ; i < N-1 ; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(pts_b.row(i),h_sq,i+3*N-2) << std::endl;
  }

  // Adding the corresponding lines for the square
  add_lines(out, 4*N-4, 1, 1);

  // Adding points corresponding to kite
  // # of points on the kite
  unsigned N_kite = M;
  //Eigen::VectorXd temp = Eigen::VectorXd::LinSpaced(N_kite + 1, 0, 2 * M_PI);
  //Eigen::VectorXd theta = temp.head(N_kite);
  Eigen::VectorXd theta = get_kite_params(N_kite);
  for (unsigned i = 0; i < N_kite; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(kite(theta(i)), h_kite, 4*N-3 + i) << std::endl;
  }
  add_lines(out, N_kite, 4*N-3, 4*N-3);

  // Adding the curve loops
  add_curve_loop(out,1,1,4*N-4);
  add_physical_curve(out,1,1,4*N-4,"Outer boundary");
  add_curve_loop(out,2,4*N-3,N_kite);
  add_physical_curve(out,2,4*N-3,N_kite,"Inner boundary");
  // Adding the plane surface
  out << new_entity_command << std::endl;
  out << "Plane Surface(1) = {1,2};" << std::endl;
  //std::cout << get_kite_params(5) << std::endl;
  return 0;
}
