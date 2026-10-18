#include "dqmc/current.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

bool finite(double value) { return std::isfinite(value); }

void check(bool ok, const std::string& message) {
  if (!ok) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}

void check_near(double actual, double expected, double tolerance,
                const std::string& message) {
  if (!finite(actual) || !finite(expected) || !finite(tolerance) ||
      tolerance < 0.0) {
    ++failures;
    std::cerr << "FAIL: " << message << ": non-finite value or invalid tolerance\n";
    return;
  }
  if (std::abs(actual - expected) > tolerance) {
    ++failures;
    std::cerr.precision(14);
    std::cerr << "FAIL: " << message << ": actual=" << actual
              << " expected=" << expected << " tolerance=" << tolerance << '\n';
  }
}

void check_matrix_near(const dqmc::Matrix& actual, const dqmc::Matrix& expected,
                       double tolerance, const std::string& message) {
  check(actual.allFinite() && expected.allFinite(), message + ": non-finite matrix");
  check_near((actual - expected).norm(), 0.0, tolerance, message);
}

void test_inverse_identity(const dqmc::Matrix& matrix, const std::string& name) {
  const dqmc::Matrix identity = dqmc::Matrix::Identity(matrix.rows(), matrix.cols());
  check_matrix_near(dqmc::test_support::inverse_identity_plus(matrix),
                    (identity + matrix).inverse(), 1e-11,
                    "inverse_identity_plus " + name);
}

void test_inverse_pair(const dqmc::Matrix& u, const dqmc::Matrix& v,
                       const std::string& name) {
  check_matrix_near(dqmc::test_support::inverse_inverse_plus(u, v),
                    (u.inverse() + v).inverse(), 1e-11,
                    "inverse_inverse_plus " + name);
}

void test_inverse_helpers() {
  const dqmc::Matrix identity = dqmc::Matrix::Identity(3, 3);
  const dqmc::Matrix diagonal =
      (dqmc::Vector(3) << 2.0, 3.0, 4.0).finished().asDiagonal();
  dqmc::Matrix permutation(3, 3);
  permutation << 0, 1, 0, 0, 0, 1, 1, 0, 0;
  dqmc::Matrix a(3, 3), b(3, 3);
  a << 2.0, 0.2, 0.1, 0.3, 3.0, 0.4, 0.1, 0.2, 4.0;
  b << 1.5, 0.1, 0.0, 0.2, 2.0, 0.3, 0.1, 0.0, 2.5;

  test_inverse_identity(identity, "identity");
  test_inverse_identity(diagonal, "diagonal");
  test_inverse_identity(permutation, "permutation");
  test_inverse_identity(a, "well-conditioned matrix");
  test_inverse_pair(identity, identity, "identity pair");
  test_inverse_pair(diagonal, diagonal, "diagonal pair");
  test_inverse_pair(permutation, diagonal, "permuted pair");
  test_inverse_pair(a, b, "well-conditioned pair");

  std::vector<dqmc::Matrix> factors;
  factors.reserve(40);
  const dqmc::Matrix factor_a =
      (dqmc::Matrix(3, 3) << 0.50, 0.08, 0.03, 0.04, 0.60, 0.06, 0.02, 0.05, 0.70)
          .finished();
  const dqmc::Matrix factor_b =
      (dqmc::Matrix(3, 3) << 0.55, 0.04, 0.02, 0.03, 0.65, 0.05, 0.01, 0.04, 0.60)
          .finished();
  dqmc::Matrix direct = dqmc::Matrix::Identity(3, 3);
  for (int i = 0; i < 40; ++i) {
    const auto& factor = i % 2 == 0 ? factor_a : factor_b;
    factors.push_back(factor);
    direct = factor * direct;
  }
  check_matrix_near(dqmc::test_support::product_ldr(factors), direct, 1e-10,
                    "40-factor LDR reconstruction");
}

void test_exhaustive_response() {
  for (int slices = 1; slices <= 4; ++slices) {
    dqmc::Parameters params;
    params.lx = 3;
    params.ly = 1;
    params.slices = slices;
    params.abs_u = 2.0;
    params.mu = -0.3;
    params.dtau = 0.2;
    dqmc::Solver solver(params);
    double partition = 0.0;
    double response = 0.0;
    const int field_count = 1 << (3 * slices);
    for (int bits = 0; bits < field_count; ++bits) {
      std::vector<std::int8_t> fields(3 * slices);
      for (int i = 0; i < 3 * slices; ++i)
        fields[i] = ((bits >> i) & 1) ? 1 : -1;
      solver.set_fields_for_test(fields);
      const double weight = solver.field_weight_for_test();
      check(finite(weight) && weight >= 0.0, "finite nonnegative field weight");
      partition += weight;
      response += weight * solver.current_response(1, 0).response;
    }
    const auto exact = dqmc::exact_tiny_current_response(
        solver.kinetic(), params.abs_u, slices * params.dtau, params.t,
        params.lx, params.ly, 1, 0);
    check(finite(partition) && partition > 0.0, "positive exhaustive partition");
    check_near(response / partition, exact.response, 0.06,
               "exhaustive response M=" + std::to_string(slices));
  }
}

void test_longitudinal_gauge() {
  dqmc::Parameters params;
  params.lx = 3;
  params.ly = 1;
  params.slices = 40;
  params.abs_u = 2.0;
  params.mu = -0.3;
  params.dtau = 0.1;
  params.seed = 881;
  const auto exact = dqmc::exact_tiny_current_response(
      dqmc::square_lattice_kinetic(params.lx, params.ly, params.t, params.mu),
      params.abs_u, params.slices * params.dtau, params.t,
      params.lx, params.ly, 1, 0);
  check_near(exact.response, exact.minus_kinetic_x, 2e-9,
             "exact nonzero-longitudinal gauge identity");

  dqmc::Solver solver(params);
  for (int i = 0; i < 2000; ++i) solver.sweep();
  double response = 0.0;
  double kinetic = 0.0;
  constexpr int count = 10000;
  for (int i = 0; i < count; ++i) {
    solver.sweep();
    const auto current = solver.current_response(1, 0);
    response += current.response;
    kinetic += current.minus_kinetic_x;
  }
  response /= count;
  kinetic /= count;
  check_near(response, kinetic, 0.03, "DQMC nonzero-longitudinal gauge identity");
  check_near(response, exact.response, 0.04,
             "DQMC nonzero-longitudinal response vs exact Kubo reference");
  check_near(kinetic, exact.minus_kinetic_x, 0.03,
             "DQMC longitudinal diamagnetic term vs exact reference");
  std::cout << "longitudinal exact=" << exact.response
            << " DQMC=" << response << " kinetic=" << kinetic << '\n';
}

void test_transverse_response() {
  dqmc::Parameters params;
  params.lx = 2;
  params.ly = 2;
  params.slices = 40;
  params.abs_u = 2.0;
  params.mu = -0.3;
  params.dtau = 0.1;
  params.seed = 918;
  const auto kinetic_matrix =
      dqmc::square_lattice_kinetic(params.lx, params.ly, params.t, params.mu);
  const auto exact = dqmc::exact_tiny_current_response(
      kinetic_matrix, params.abs_u, params.slices * params.dtau, params.t,
      params.lx, params.ly, 0, 1);
  const auto exact_nonzero = dqmc::exact_tiny_current_response(
      kinetic_matrix, params.abs_u, params.slices * params.dtau, params.t,
      params.lx, params.ly, 1, 0);
  check(std::abs(exact_nonzero.response) > 1e-6,
        "2D nonzero-momentum exact reference is nonzero");

  dqmc::Solver solver(params);
  for (int i = 0; i < 2000; ++i) solver.sweep();
  double response = 0.0;
  double response_nonzero = 0.0;
  double kinetic = 0.0;
  constexpr int count = 10000;
  for (int i = 0; i < count; ++i) {
    solver.sweep();
    const auto current = solver.current_response(0, 1);
    response += current.response;
    response_nonzero += solver.current_response(1, 0).response;
    kinetic += current.minus_kinetic_x;
  }
  response /= count;
  response_nonzero /= count;
  kinetic /= count;
  check_near(response, exact.response, 0.06,
             "DQMC 2D transverse response vs exact Kubo reference");
  check_near(response_nonzero, exact_nonzero.response, 0.06,
             "DQMC nonzero 2D momentum response vs exact Kubo reference");
  check_near(kinetic, exact.minus_kinetic_x, 0.03,
             "DQMC 2D diamagnetic term vs exact reference");
  std::cout << "2D transverse exact=" << exact.response
            << " DQMC=" << response << "; 2D nonzero exact="
            << exact_nonzero.response << " DQMC=" << response_nonzero
            << "; kinetic=" << kinetic << '\n';
}

}  // namespace

int main() {
  test_inverse_helpers();
  test_exhaustive_response();
  test_longitudinal_gauge();
  test_transverse_response();
  if (failures == 0) std::cout << "All current-response tests passed.\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
