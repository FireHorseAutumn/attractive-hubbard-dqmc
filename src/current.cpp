#include "dqmc/current.hpp"
#include "dqmc/ldr.hpp"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

namespace dqmc {
namespace {

using Complex = std::complex<double>;
using ComplexMatrix = Eigen::MatrixXcd;

ComplexMatrix current_matrix(int lx, int ly, double t, int qx, int qy) {
  const int n = lx * ly;
  ComplexMatrix j = ComplexMatrix::Zero(n, n);
  const Complex imaginary(0.0, 1.0);
  const double pi = std::acos(-1.0);
  for (int y = 0; y < ly; ++y) {
    for (int x = 0; x < lx; ++x) {
      const int site = y * lx + x;
      const int right = y * lx + (x + 1) % lx;
      const double phase_real = 2.0 * pi *
          (static_cast<double>(qx * x) / lx + static_cast<double>(qy * y) / ly);
      const Complex phase = std::exp(imaginary * phase_real);
      j(right, site) += imaginary * t * phase;
      j(site, right) -= imaginary * t * phase;
    }
  }
  return j;
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
  if (parity(state & (source_bit - 1))) sign = -sign;
  const int next = state ^ source_bit;
  if (parity(next & (destination_bit - 1))) sign = -sign;
  return {next | destination_bit, static_cast<double>(sign)};
}

void add_current_to_fock(ComplexMatrix& j, int state, int nsites, int lx,
                         double t, int qx, int qy) {
  const int ly = nsites / lx;
  const Complex imaginary(0.0, 1.0);
  const double pi = std::acos(-1.0);
  for (int y = 0; y < ly; ++y) {
    for (int x = 0; x < lx; ++x) {
      const int site = y * lx + x;
      const int right = y * lx + (x + 1) % lx;
      const double angle = 2.0 * pi *
          (static_cast<double>(qx * x) / lx + static_cast<double>(qy * y) / ly);
      const Complex phase = std::exp(imaginary * angle);
      for (int spin = 0; spin < 2; ++spin) {
        const HopResult forward = hop(state, 2 * right + spin, 2 * site + spin);
        if (forward.sign != 0.0)
          j(forward.state, state) += imaginary * t * phase * forward.sign;
        const HopResult backward = hop(state, 2 * site + spin, 2 * right + spin);
        if (backward.sign != 0.0)
          j(backward.state, state) -= imaginary * t * phase * backward.sign;
      }
    }
  }
}

}  // namespace

namespace test_support {

Matrix inverse_identity_plus(const Matrix& matrix) {
  return detail::inverse_identity_plus(detail::make_ldr(matrix));
}

Matrix inverse_inverse_plus(const Matrix& u, const Matrix& v) {
  return detail::inverse_inverse_plus(detail::make_ldr(u), detail::make_ldr(v));
}

Matrix product_ldr(const std::vector<Matrix>& factors) {
  if (factors.empty()) throw std::invalid_argument("empty LDR product");
  auto result = detail::identity_ldr(factors.front().rows());
  for (const auto& factor : factors) result = detail::prepend_ldr(result, factor);
  return result.left * result.diagonal.asDiagonal() * result.right;
}

}  // namespace test_support

double Solver::minus_kinetic_x() const {
  const Matrix equal = green_at_slice(0);
  const Matrix occupations = Matrix::Identity(nsites_, nsites_) - equal.transpose();
  double value = 0.0;
  for (int y = 0; y < params_.ly; ++y) {
    for (int x = 0; x < params_.lx; ++x) {
      const int site = y * params_.lx + x;
      const int right = y * params_.lx + (x + 1) % params_.lx;
      value += 2.0 * params_.t *
               (occupations(right, site) + occupations(site, right));
    }
  }
  return value / nsites_;
}

CurrentResult Solver::current_response(int qx, int qy) const {
  if (qx < 0 || qx >= params_.lx || qy < 0 || qy >= params_.ly)
    throw std::out_of_range("momentum index outside lattice");
  const auto identity = detail::identity_ldr(nsites_);
  std::vector<detail::Ldr> p_ldr(params_.slices + 1, identity);
  std::vector<detail::Ldr> q_ldr(params_.slices + 1, identity);
  for (int l = 0; l < params_.slices; ++l)
    p_ldr[l + 1] = detail::prepend_ldr(p_ldr[l], bmat_[l]);
  for (int l = 0; l < params_.slices; ++l) {
    for (int r = l; r < params_.slices; ++r)
      q_ldr[l] = detail::prepend_ldr(q_ldr[l], bmat_[r]);
  }
  const ComplexMatrix j = current_matrix(params_.lx, params_.ly, params_.t, qx, qy);
  const ComplexMatrix j_adjoint = j.adjoint();
  double response = 0.0;
  for (int l = 0; l < params_.slices; ++l) {
    const Matrix g_tau0 = detail::inverse_inverse_plus(p_ldr[l], q_ldr[l]);
    const Matrix g_zero_minus = -detail::inverse_inverse_plus(q_ldr[l], p_ldr[l]);
    const Complex connected =
        -(j * g_tau0 * j_adjoint * g_zero_minus).trace();
    response += 2.0 * params_.dtau * connected.real() / nsites_;
  }
  return {response, minus_kinetic_x()};
}

CurrentResult exact_tiny_current_response(const Matrix& kinetic, double abs_u,
                                          double beta, double t, int lx, int ly,
                                          int qx, int qy) {
  const int nsites = static_cast<int>(kinetic.rows());
  if (nsites > 4 || nsites != lx * ly)
    throw std::invalid_argument("exact current reference requires <=4 sites");
  const int norb = 2 * nsites;
  const int dimension = 1 << norb;
  Matrix hamiltonian = Matrix::Zero(dimension, dimension);
  ComplexMatrix j = ComplexMatrix::Zero(dimension, dimension);
  Matrix minus_kx = Matrix::Zero(dimension, dimension);

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
        for (int k = 0; k < nsites; ++k) {
          if (i == k || kinetic(i, k) == 0.0) continue;
          const HopResult h = hop(state, 2 * i + spin, 2 * k + spin);
          if (h.sign != 0.0)
            hamiltonian(h.state, state) += kinetic(i, k) * h.sign;
        }
      }
    }
    add_current_to_fock(j, state, nsites, lx, t, qx, qy);

    for (int y = 0; y < ly; ++y) {
      for (int x = 0; x < lx; ++x) {
        const int site = y * lx + x;
        const int right = y * lx + (x + 1) % lx;
        for (int spin = 0; spin < 2; ++spin) {
          const HopResult forward = hop(state, 2 * right + spin, 2 * site + spin);
          if (forward.sign != 0.0)
            minus_kx(forward.state, state) += t * forward.sign;
          const HopResult backward = hop(state, 2 * site + spin, 2 * right + spin);
          if (backward.sign != 0.0)
            minus_kx(backward.state, state) += t * backward.sign;
        }
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(hamiltonian);
  if (eigensolver.info() != Eigen::Success)
    throw std::runtime_error("exact current diagonalization failed");
  const Vector energies = eigensolver.eigenvalues();
  const Matrix vectors = eigensolver.eigenvectors();
  const Vector weights = (-beta * energies.array()).exp();
  const Vector probabilities = weights / weights.sum();
  const ComplexMatrix j_eigen = vectors.adjoint() * j * vectors;
  double response = 0.0;
  for (int n = 0; n < dimension; ++n) {
    for (int m = 0; m < dimension; ++m) {
      const double gap = energies(m) - energies(n);
      double factor = 0.0;
      if (std::abs(gap) < 1e-10) {
        factor = beta * probabilities(n);
      } else {
        factor = (probabilities(n) - probabilities(m)) / gap;
      }
      response += std::norm(j_eigen(n, m)) * factor / nsites;
    }
  }
  const Matrix density_matrix =
      vectors * probabilities.asDiagonal() * vectors.transpose();
  const double kinetic_x = (density_matrix * minus_kx).trace() / nsites;
  return {response, kinetic_x};
}

}  // namespace dqmc
