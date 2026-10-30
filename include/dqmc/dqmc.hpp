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

  // Checkpoints are valid at sweep boundaries. Saving requires at least one
  // recompute flag so the boundary Green function is a pure function of the
  // stored state. The format is native-endian and versioned; a trailing
  // checksum covers the whole payload.
  static constexpr std::uint32_t checkpoint_version = 1;
  void save_checkpoint(const std::string& path, std::uint64_t sweep_tag) const;
  static Solver load_checkpoint(const std::string& path);
  std::uint64_t sweep_tag() const { return sweep_tag_; }

  Matrix green_at_slice(int slice) const;
  double determinant_at_slice_zero() const;
  void set_fields_for_test(const std::vector<std::int8_t>& fields);
  double field_weight_for_test() const;
  const Matrix& green() const { return green_; }
  const Matrix& kinetic() const { return kinetic_; }
  double lambda() const { return lambda_; }
  std::uint64_t accepted() const { return accepted_; }
  std::uint64_t attempted() const { return attempted_; }
  const Parameters& params() const { return params_; }

private:
  struct UninitializedTag {};
  Solver(Parameters params, UninitializedTag);
  void build_derived();

  Parameters params_;
  int nsites_;
  double lambda_;
  // Ising fields are ±1, so every field-factor exponent is -2λ, 0, or +2λ;
  // cache the three values once (same std::exp calls, same arguments as the
  // per-element form — bitwise identical).
  double exp_m2lambda_ = 0.0;
  double exp_0lambda_ = 1.0;
  double exp_p2lambda_ = 0.0;
  Matrix kinetic_;
  Matrix exp_k_;
  Matrix exp_minus_k_;
  Matrix green_;
  Matrix scaled_scratch_;  // reusable per-slice rescale target
  Vector rk_col_;          // reusable rank-1 update scratch
  Vector rk_row_;
  std::vector<Matrix> bmat_;
  std::vector<std::int8_t> fields_;
  std::mt19937_64 rng_;
  std::uint64_t accepted_ = 0;
  std::uint64_t attempted_ = 0;
  std::uint64_t sweep_tag_ = 0;

  void rebuild_bmatrices();
  void propagate_to_next_slice(int slice);
  double attempt_flip(int site, int slice);
};

}  // namespace dqmc
