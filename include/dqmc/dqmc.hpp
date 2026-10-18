#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <random>
#include <vector>

namespace dqmc {

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

struct ReferenceResult {
  double density = 0.0;
  double energy = 0.0;
  double double_occupancy = 0.0;
  double pair_structure = 0.0;
  Matrix green;
};

struct CalibrationResult {
  double mu = 0.0;
  double density = 0.0;
  int iterations = 0;
};

struct CurrentResult {
  double response = 0.0;
  double minus_kinetic_x = 0.0;
};

struct Parameters {
  int lx = 2;
  int ly = 1;
  int slices = 20;
  double t = 1.0;
  double abs_u = 2.0;
  double mu = 0.0;
  double dtau = 0.1;
  std::uint64_t seed = 12345;
  bool recompute_each_slice = true;
  bool recompute_after_sweep = true;
};

Matrix square_lattice_kinetic(int lx, int ly, double t, double mu);
double hs_lambda(double dtau, double abs_u);
ReferenceResult free_fermion_reference(const Matrix& kinetic, double beta);
ReferenceResult exact_tiny_reference(const Matrix& kinetic, double abs_u,
                                     double beta);
CalibrationResult calibrate_mu(Parameters params, double target_density,
                               int warmup, int measure, double tolerance,
                               int max_iterations);

class Solver {
public:
  explicit Solver(Parameters params);

  void initialize();
  void sweep();
  double flip_for_test(int site);
  double density() const;
  double energy() const;
  double double_occupancy() const;
  double pair_structure() const;
  ReferenceResult measure() const;
  double minus_kinetic_x() const;
  CurrentResult current_response(int qx, int qy) const;

  Matrix green_at_slice(int slice) const;
  double determinant_at_slice_zero() const;
  void set_fields_for_test(const std::vector<std::int8_t>& fields);
  double field_weight_for_test() const;
  const Matrix& green() const { return green_; }
  const Matrix& kinetic() const { return kinetic_; }
  double lambda() const { return lambda_; }
  std::uint64_t accepted() const { return accepted_; }
  std::uint64_t attempted() const { return attempted_; }

private:
  Parameters params_;
  int nsites_;
  double lambda_;
  Matrix kinetic_;
  Matrix exp_k_;
  Matrix exp_minus_k_;
  Matrix green_;
  std::vector<Matrix> bmat_;
  std::vector<std::int8_t> fields_;
  std::mt19937_64 rng_;
  std::uint64_t accepted_ = 0;
  std::uint64_t attempted_ = 0;

  void rebuild_bmatrices();
  void propagate_to_next_slice(int slice);
  double attempt_flip(int site, int slice);
};

}  // namespace dqmc
