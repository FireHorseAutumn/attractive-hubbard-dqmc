#include "dqmc/dqmc.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

int failures = 0;

void check(bool condition, const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}

void check_near(double actual, double expected, double tolerance,
                const std::string& message) {
  if (std::abs(actual - expected) > tolerance) {
    ++failures;
    std::cerr.precision(14);
    std::cerr << "FAIL: " << message << ": actual=" << actual
              << " expected=" << expected << " tolerance=" << tolerance << '\n';
  }
}

void test_hs_identity() {
  for (int up = 0; up <= 1; ++up) {
    for (int down = 0; down <= 1; ++down) {
      for (double a : {0.0, 0.05, 0.2, 0.7}) {
        const double lambda = std::acosh(std::exp(a / 2.0));
        const double lhs = std::exp(a * (up - 0.5) * (down - 0.5));
        double rhs = 0.0;
        for (double s : {-1.0, 1.0})
          rhs += 0.5 * std::exp(-a / 4.0) *
                 std::exp(lambda * s * (up + down - 1.0));
        check_near(rhs, lhs, 2e-14, "four-state HS identity");
      }
    }
  }
}

void test_free_fermion() {
  dqmc::Parameters params;
  params.lx = 2;
  params.ly = 2;
  params.slices = 8;
  params.abs_u = 0.0;
  params.mu = -0.3;
  params.dtau = 0.1;
  dqmc::Solver solver(params);
  const auto exact = dqmc::free_fermion_reference(solver.kinetic(),
                                                  params.slices * params.dtau);
  const auto measured = solver.measure();
  check((solver.green() - exact.green).norm() < 1e-11,
        "free-fermion Green function");
  check_near(measured.density, exact.density, 1e-11, "free-fermion density");
  check_near(measured.energy, exact.energy, 1e-11, "free-fermion energy");
  check_near(measured.double_occupancy, exact.double_occupancy, 1e-11,
             "free-fermion double occupancy");
}

void test_local_update() {
  dqmc::Parameters params;
  params.lx = 2;
  params.ly = 2;
  params.slices = 5;
  params.abs_u = 2.0;
  params.mu = -0.2;
  params.dtau = 0.15;
  params.recompute_after_sweep = false;
  dqmc::Solver solver(params);

  const double old_det = solver.determinant_at_slice_zero();
  const double ratio = solver.flip_for_test(1);
  const double new_det = solver.determinant_at_slice_zero();
  check_near(new_det / old_det, ratio, 2e-10,
             "determinant ratio for forced local flip");
  const dqmc::Matrix full = solver.green_at_slice(0);
  check((solver.green() - full).norm() < 2e-10,
        "Green update for forced local flip");

  solver.sweep();
  const dqmc::Matrix after_sweeps = solver.green_at_slice(0);
  check((solver.green() - after_sweeps).norm() < 2e-8,
        "cyclic Green propagation across one unstabilized sweep");
}

void test_tiny_interacting() {
  dqmc::Parameters params;
  params.lx = 2;
  params.ly = 1;
  params.slices = 40;
  params.abs_u = 2.0;
  params.mu = -0.3;
  params.dtau = 0.1;
  params.seed = 777;
  dqmc::Solver solver(params);
  const auto exact = dqmc::exact_tiny_reference(
      solver.kinetic(), params.abs_u, params.slices * params.dtau);

  for (int i = 0; i < 2000; ++i) solver.sweep();
  dqmc::ReferenceResult average;
  constexpr int count = 40000;
  for (int i = 0; i < count; ++i) {
    solver.sweep();
    const auto current = solver.measure();
    average.density += current.density;
    average.energy += current.energy;
    average.double_occupancy += current.double_occupancy;
    average.pair_structure += current.pair_structure;
  }
  average.density /= count;
  average.energy /= count;
  average.double_occupancy /= count;
  average.pair_structure /= count;

  check_near(average.density, exact.density, 0.01,
             "tiny interacting DQMC density vs exact diagonalization");
  check_near(average.energy, exact.energy, 0.04,
             "tiny interacting DQMC energy vs exact diagonalization");
  check_near(average.double_occupancy, exact.double_occupancy, 0.01,
             "tiny interacting DQMC double occupancy vs exact diagonalization");
  check_near(average.pair_structure, exact.pair_structure, 0.02,
             "tiny interacting DQMC pair structure vs exact diagonalization");
  std::cout << "tiny exact: density=" << exact.density
            << " energy=" << exact.energy
            << " double=" << exact.double_occupancy
            << " pair=" << exact.pair_structure << '\n'
            << "tiny DQMC: density=" << average.density
            << " energy=" << average.energy
            << " double=" << average.double_occupancy
            << " pair=" << average.pair_structure << '\n';
}

void test_calibration() {
  dqmc::Parameters params;
  params.lx = 2;
  params.ly = 1;
  params.slices = 20;
  params.abs_u = 1.0;
  params.dtau = 0.1;
  params.seed = 9182;
  constexpr double target = 1.2;
  const auto result = dqmc::calibrate_mu(params, target, 100, 300, 0.03, 14);
  check(std::abs(result.density - target) <= 0.03,
        "chemical-potential calibration return value");

  params.mu = result.mu;
  params.seed += 1;
  dqmc::Solver solver(params);
  for (int i = 0; i < 1000; ++i) solver.sweep();
  double density = 0.0;
  for (int i = 0; i < 10000; ++i) {
    solver.sweep();
    density += solver.density();
  }
  density /= 10000;
  check_near(density, target, 0.08,
             "independent run at calibrated chemical potential");
  std::cout << "calibration: mu=" << result.mu << " density=" << density
            << " iterations=" << result.iterations << '\n';
}

}  // namespace

int main() {
  test_hs_identity();
  test_free_fermion();
  test_local_update();
  test_tiny_interacting();
  test_calibration();
  if (failures == 0) std::cout << "All core DQMC tests passed.\n";
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
