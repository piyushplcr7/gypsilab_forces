#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

Eigen::Vector2d kite(double t) {
  return Eigen::Vector2d(.35 * std::cos(t) + .1625 * std::cos(2 * t),
                         .35 * std::sin(t));
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

int main() {
  // Creating the geometry file
  std::string filename = "kite_sq.geo";
  std::ofstream out(filename);
  std::string new_entity_command = "//+";

  // Mesh size parameters
  double h_sq = 1;
  double h_kite = 1;

  // Edge length of the square
  double s = 6;

  // Adding points of the square
  out << new_entity_command << std::endl;
  out << add_point(Eigen::Vector2d(3, 3), h_sq, 1) << std::endl;
  out << new_entity_command << std::endl;
  out << add_point(Eigen::Vector2d(-3, 3), h_sq, 2) << std::endl;
  out << new_entity_command << std::endl;
  out << add_point(Eigen::Vector2d(-3, -3), h_sq, 3) << std::endl;
  out << new_entity_command << std::endl;
  out << add_point(Eigen::Vector2d(3, -3), h_sq, 4) << std::endl;

  // Adding the corresponding lines
  add_lines(out, 4, 1, 1);

  // Adding points corresponding to kite
  // # of points on the kite
  unsigned N = 5;
  Eigen::VectorXd temp = Eigen::VectorXd::LinSpaced(N + 1, 0, 2 * M_PI);
  Eigen::VectorXd theta = temp.head(N);
  for (unsigned i = 0; i < N; ++i) {
    out << new_entity_command << std::endl;
    out << add_point(kite(theta(i)), h_sq, 5 + i) << std::endl;
  }
  add_lines(out, N, 5, 5);

  return 0;
}
