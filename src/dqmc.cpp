#include "dqmc/dqmc.hpp"
#include "dqmc/ldr.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace dqmc {
namespace {

int orbital_count(int state_count) {
  int n = 0;
  while ((1 << n) < state_count) ++n;
  return n;
}

struct HopResult {
  int state;
  double sign;
};

HopResult hop(int state, int destination, int source) {
  const int source_bit = 1 << source;
  const int destination_bit = 1 << destination;
  if ((state & source_bit) == 0 || (state & destination_bit) != 0)
    return {state, 0.0};
  const auto parity = [](int value) { return __builtin_parity(value); };
  int sign = 1;
  const int source_below = state & (source_bit - 1);
  if (parity(source_below)) sign = -sign;
  int next = state ^ source_bit;
  const int destination_below = next & (destination_bit - 1);
  if (parity(destination_below)) sign = -sign;
  next |= destination_bit;
  return {next, static_cast<double>(sign)};
}

HopResult annihilate(int state, int orbital) {
  const int bit = 1 << orbital;
  if ((state & bit) == 0) return {state, 0.0};
  const auto parity = [](int value) { return __builtin_parity(value); };
  const double sign = parity(state & (bit - 1)) ? -1.0 : 1.0;
  return {state ^ bit, sign};
}

HopResult annihilate_pair(int state, int site) {
  const HopResult up = annihilate(state, 2 * site);
  if (up.sign == 0.0) return {state, 0.0};
  const HopResult down = annihilate(up.state, 2 * site + 1);
  return {down.state, up.sign * down.sign};
}

}  // namespace

Matrix square_lattice_kinetic(int lx, int ly, double t, double mu) {
  if (lx <= 0 || ly <= 0) throw std::invalid_argument("invalid lattice size");
  const int n = lx * ly;
  Matrix k = Matrix::Zero(n, n);
  for (int y = 0; y < ly; ++y) {
    for (int x = 0; x < lx; ++x) {
      const int i = y * lx + x;
      const int right = y * lx + (x + 1) % lx;
      const int down = ((y + 1) % ly) * lx + x;
      if (right != i) {
        k(i, right) -= t;
        k(right, i) -= t;
      }
      if (down != i) {
        k(i, down) -= t;
        k(down, i) -= t;
      }
    }
  }
  k.diagonal().array() -= mu;
  return k;
}

double hs_lambda(double dtau, double abs_u) {
  return std::acosh(std::exp(0.5 * dtau * abs_u));
}

ReferenceResult free_fermion_reference(const Matrix& kinetic, double beta) {
  Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(kinetic);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error("free-fermion diagonalization failed");
  const Vector energies = eigensolver.eigenvalues();
  const Matrix vectors = eigensolver.eigenvectors();
  ReferenceResult result;
  Vector green_diagonal(energies.size());
  for (int i = 0; i < energies.size(); ++i) {
    const double f = 1.0 / (1.0 + std::exp(beta * energies(i)));
    green_diagonal(i) = 1.0 - f;
    result.density += 2.0 * f;
    result.energy += 2.0 * energies(i) * f;
  }
  result.density /= kinetic.rows();
  result.green = vectors * green_diagonal.asDiagonal() * vectors.transpose();
  double per_spin_density = 0.0;
  for (int i = 0; i < kinetic.rows(); ++i)
    per_spin_density += 1.0 - result.green(i, i);
  per_spin_density /= kinetic.rows();
  result.double_occupancy = per_spin_density * per_spin_density;
  return result;
}

ReferenceResult exact_tiny_reference(const Matrix& kinetic, double abs_u,
                                     double beta) {
  const int nsites = static_cast<int>(kinetic.rows());
  if (nsites > 4) throw std::invalid_argument("exact reference limited to 4 sites");
  const int norb = 2 * nsites;
  const int dimension = 1 << norb;
  if (orbital_count(dimension) != norb)
    throw std::logic_error("invalid exact-reference dimension");

  Matrix hamiltonian = Matrix::Zero(dimension, dimension);
  for (int state = 0; state < dimension; ++state) {
    double diagonal = 0.0;
    for (int site = 0; site < nsites; ++site) {
      const int up = 2 * site;
      const int down = up + 1;
      const double nup = (state & (1 << up)) ? 1.0 : 0.0;
      const double ndown = (state & (1 << down)) ? 1.0 : 0.0;
      diagonal += kinetic(site, site) * (nup + ndown);
      diagonal += -abs_u * (nup - 0.5) * (ndown - 0.5);
    }
    hamiltonian(state, state) += diagonal;

    for (int spin = 0; spin < 2; ++spin) {
      for (int i = 0; i < nsites; ++i) {
        for (int j = 0; j < nsites; ++j) {
          if (i == j || kinetic(i, j) == 0.0) continue;
          const int destination = 2 * i + spin;
          const int source = 2 * j + spin;
          const HopResult h = hop(state, destination, source);
          if (h.sign != 0.0)
            hamiltonian(h.state, state) += kinetic(i, j) * h.sign;
        }
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(hamiltonian);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error("tiny exact diagonalization failed");
  const Vector energies = eigensolver.eigenvalues();
  const Matrix vectors = eigensolver.eigenvectors();
  const Vector weights = (-beta * energies.array()).exp();
  const double partition = weights.sum();

  const Matrix density_matrix =
      vectors * (weights / partition).asDiagonal() * vectors.transpose();
  Matrix pair_sum = Matrix::Zero(dimension, dimension);
  ReferenceResult result;
  for (int state = 0; state < dimension; ++state) {
    double number = 0.0;
    double double_sum = 0.0;
    for (int site = 0; site < nsites; ++site) {
      const bool up = (state & (1 << (2 * site))) != 0;
      const bool down = (state & (1 << (2 * site + 1))) != 0;
      number += static_cast<double>(up) + static_cast<double>(down);
      double_sum += (up && down) ? 1.0 : 0.0;
      const HopResult pair = annihilate_pair(state, site);
      if (pair.sign != 0.0)
        pair_sum(pair.state, state) += pair.sign / std::sqrt(nsites);
    }
    const double eigen_weight =
        vectors.row(state).array().square().matrix().dot(weights) / partition;
    result.density += eigen_weight * number / nsites;
    result.double_occupancy += eigen_weight * double_sum / nsites;
  }
  result.energy = energies.dot(weights) / partition;
  const Matrix pair_fluctuation = pair_sum.adjoint() * pair_sum +
                                  pair_sum * pair_sum.adjoint();
  result.pair_structure = (density_matrix * pair_fluctuation).trace();
  return result;
}

Solver::Solver(Parameters params)
    : params_(params),
      nsites_(params.lx * params.ly),
      lambda_(hs_lambda(params.dtau, params.abs_u)),
      kinetic_(square_lattice_kinetic(params.lx, params.ly, params.t, params.mu)),
      rng_(params.seed) {
  build_derived();
  initialize();
}

Solver::Solver(Parameters params, UninitializedTag)
    : params_(params),
      nsites_(params.lx * params.ly),
      lambda_(hs_lambda(params.dtau, params.abs_u)),
      kinetic_(square_lattice_kinetic(params.lx, params.ly, params.t, params.mu)),
      rng_(params.seed) {
  build_derived();
}

void Solver::build_derived() {
  if (nsites_ <= 0 || params_.slices <= 0 || params_.dtau <= 0.0)
    throw std::invalid_argument("invalid DQMC parameters");
  Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(kinetic_);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error("kinetic diagonalization failed");
  const Vector values = eigensolver.eigenvalues();
  const Matrix vectors = eigensolver.eigenvectors();
  exp_k_ = vectors * (params_.dtau * values.array()).exp().matrix().asDiagonal() *
           vectors.transpose();
  exp_minus_k_ = vectors *
                 (-params_.dtau * values.array()).exp().matrix().asDiagonal() *
                 vectors.transpose();
}

void Solver::initialize() {
  fields_.resize(nsites_ * params_.slices);
  std::uniform_int_distribution<int> sign(0, 1);
  for (auto& field : fields_) field = sign(rng_) == 0 ? -1 : 1;
  rebuild_bmatrices();
  green_ = green_at_slice(0);
  accepted_ = 0;
  attempted_ = 0;
}

void Solver::rebuild_bmatrices() {
  bmat_.assign(params_.slices, Matrix::Zero(nsites_, nsites_));
  for (int l = 0; l < params_.slices; ++l) {
    bmat_[l] = exp_minus_k_;
    for (int i = 0; i < nsites_; ++i)
      bmat_[l].col(i) *= std::exp(lambda_ * fields_[l * nsites_ + i]);
  }
}

Matrix Solver::green_at_slice(int slice) const {
  if (slice < 0 || slice >= params_.slices)
    throw std::out_of_range("invalid slice");
  auto product = detail::identity_ldr(nsites_);
  for (int l = slice; l < params_.slices; ++l)
    product = detail::prepend_ldr(product, bmat_[l]);
  for (int l = 0; l < slice; ++l)
    product = detail::prepend_ldr(product, bmat_[l]);
  return detail::inverse_identity_plus(product);
}

double Solver::determinant_at_slice_zero() const {
  Matrix product = Matrix::Identity(nsites_, nsites_);
  for (int l = 0; l < params_.slices; ++l) product = bmat_[l] * product;
  return (Matrix::Identity(nsites_, nsites_) + product).determinant();
}

void Solver::set_fields_for_test(const std::vector<std::int8_t>& fields) {
  if (static_cast<int>(fields.size()) != nsites_ * params_.slices)
    throw std::invalid_argument("field vector has wrong size");
  for (const auto field : fields) {
    if (field != -1 && field != 1)
      throw std::invalid_argument("Ising fields must be -1 or +1");
  }
  fields_ = fields;
  rebuild_bmatrices();
  green_ = green_at_slice(0);
}

double Solver::field_weight_for_test() const {
  long field_sum = 0;
  for (const auto field : fields_) field_sum += field;
  const double determinant = determinant_at_slice_zero();
  return determinant * determinant * std::exp(-lambda_ * field_sum);
}

double Solver::attempt_flip(int site, int slice) {
  const int index = slice * nsites_ + site;
  const double s = static_cast<double>(fields_[index]);
  const double r = std::exp(-2.0 * lambda_ * s);
  const double delta = r - 1.0;
  const double ratio = 1.0 + delta * (1.0 - green_(site, site));
  ++attempted_;
  const double weight_ratio = ratio * ratio * std::exp(2.0 * lambda_ * s);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  if (uniform(rng_) >= weight_ratio) return ratio;

  const Matrix update = (Matrix::Identity(nsites_, nsites_) - green_).col(site) *
                        green_.row(site);
  green_ -= (delta / ratio) * update;
  fields_[index] = static_cast<std::int8_t>(-fields_[index]);
  bmat_[slice].col(site) *= r;
  ++accepted_;
  return ratio;
}

double Solver::flip_for_test(int site) {
  const double s = static_cast<double>(fields_[site]);
  const double r = std::exp(-2.0 * lambda_ * s);
  const double delta = r - 1.0;
  const double ratio = 1.0 + delta * (1.0 - green_(site, site));
  const Matrix update = (Matrix::Identity(nsites_, nsites_) - green_).col(site) *
                        green_.row(site);
  green_ -= (delta / ratio) * update;
  fields_[site] = static_cast<std::int8_t>(-fields_[site]);
  bmat_[0].col(site) *= r;
  return ratio;
}

void Solver::propagate_to_next_slice(int slice) {
  Matrix scaled = green_;
  for (int i = 0; i < nsites_; ++i) {
    for (int j = 0; j < nsites_; ++j) {
      scaled(i, j) *= std::exp(lambda_ *
          (fields_[slice * nsites_ + i] - fields_[slice * nsites_ + j]));
    }
  }
  green_ = exp_minus_k_ * scaled * exp_k_;
}

void Solver::sweep() {
  for (int l = 0; l < params_.slices; ++l) {
    for (int i = 0; i < nsites_; ++i) attempt_flip(i, l);
    propagate_to_next_slice(l);
    if (params_.recompute_each_slice) green_ = green_at_slice((l + 1) % params_.slices);
  }
  if (params_.recompute_after_sweep) green_ = green_at_slice(0);
}

double Solver::density() const {
  double value = 0.0;
  for (int i = 0; i < nsites_; ++i) value += 2.0 * (1.0 - green_(i, i));
  return value / nsites_;
}

double Solver::double_occupancy() const {
  double value = 0.0;
  for (int i = 0; i < nsites_; ++i)
    value += (1.0 - green_(i, i)) * (1.0 - green_(i, i));
  return value / nsites_;
}

double Solver::energy() const {
  double value = 0.0;
  for (int i = 0; i < nsites_; ++i) {
    for (int j = 0; j < nsites_; ++j) {
      const double delta = i == j ? 1.0 : 0.0;
      value += 2.0 * kinetic_(i, j) * (delta - green_(j, i));
    }
  }
  for (int i = 0; i < nsites_; ++i) {
    const double nper_spin = 1.0 - green_(i, i);
    const double double_i = nper_spin * nper_spin;
    value += -params_.abs_u * (double_i - nper_spin + 0.25);
  }
  return value;
}

double Solver::pair_structure() const {
  double value = 0.0;
  for (int i = 0; i < nsites_; ++i) {
    for (int j = 0; j < nsites_; ++j) {
      const double delta = i == j ? 1.0 : 0.0;
      const double creation_annihilation = delta - green_(j, i);
      value += creation_annihilation * creation_annihilation +
               green_(i, j) * green_(i, j);
    }
  }
  return value / nsites_;
}

ReferenceResult Solver::measure() const {
  ReferenceResult result;
  result.density = density();
  result.energy = energy();
  result.double_occupancy = double_occupancy();
  result.pair_structure = pair_structure();
  result.green = green_;
  return result;
}

CalibrationResult calibrate_mu(Parameters params, double target_density,
                               int warmup, int measure, double tolerance,
                               int max_iterations) {
  if (warmup < 0 || measure <= 0 || tolerance <= 0.0 || max_iterations <= 0)
    throw std::invalid_argument("invalid calibration controls");
  const auto evaluate = [&](double mu) {
    params.mu = mu;
    Solver solver(params);
    for (int i = 0; i < warmup; ++i) solver.sweep();
    double density = 0.0;
    for (int i = 0; i < measure; ++i) {
      solver.sweep();
      density += solver.density();
    }
    return density / measure;
  };

  double low = -8.0;
  double high = 8.0;
  const double low_density = evaluate(low);
  const double high_density = evaluate(high);
  if (target_density < low_density || target_density > high_density)
    throw std::runtime_error("target density is outside the mu bracket");

  CalibrationResult result;
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const double mid = 0.5 * (low + high);
    const double density = evaluate(mid);
    result.mu = mid;
    result.density = density;
    result.iterations = iteration + 1;
    if (std::abs(density - target_density) <= tolerance) return result;
    if (density < target_density)
      low = mid;
    else
      high = mid;
  }
  throw std::runtime_error("chemical-potential calibration did not converge");
}

}  // namespace dqmc
