#include "dqmc/dqmc.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const std::string& name) {
  if (condition) {
    std::cout << "PASS " << name << "\n";
  } else {
    std::cout << "FAIL " << name << "\n";
    ++failures;
  }
}

dqmc::Parameters test_params() {
  dqmc::Parameters params;
  params.lx = 4;
  params.ly = 4;
  params.slices = 24;
  params.dtau = 0.1;
  params.abs_u = 2.0;
  params.mu = 0.3;
  params.seed = 999;
  return params;
}

struct Sample {
  double density;
  double energy;
  double double_occupancy;
  double pair_structure;
  double minus_kinetic_x;
};

Sample sample(const dqmc::Solver& solver) {
  return {solver.density(), solver.energy(), solver.double_occupancy(),
          solver.pair_structure(), solver.minus_kinetic_x()};
}

bool same_sample(const Sample& a, const Sample& b) {
  return a.density == b.density && a.energy == b.energy &&
         a.double_occupancy == b.double_occupancy &&
         a.pair_structure == b.pair_structure &&
         a.minus_kinetic_x == b.minus_kinetic_x;
}

void bitwise_trajectory(const std::string& path) {
  const int total = 30;
  const int split = 12;

  dqmc::Solver uninterrupted(test_params());
  std::vector<Sample> reference;
  for (int i = 0; i < total; ++i) {
    uninterrupted.sweep();
    reference.push_back(sample(uninterrupted));
  }

  dqmc::Solver first(test_params());
  for (int i = 0; i < split; ++i) first.sweep();
  first.save_checkpoint(path, split);

  dqmc::Solver resumed = dqmc::Solver::load_checkpoint(path);
  check(resumed.accepted() == first.accepted() &&
            resumed.attempted() == first.attempted(),
        "counters survive the round trip");
  check(resumed.sweep_tag() == static_cast<std::uint64_t>(split),
        "sweep tag survives the round trip");

  bool identical = true;
  for (int i = split; i < total; ++i) {
    resumed.sweep();
    if (!same_sample(sample(resumed), reference[static_cast<size_t>(i)]))
      identical = false;
  }
  check(identical, "resumed trajectory is bitwise identical");

  check(resumed.accepted() == uninterrupted.accepted() &&
            resumed.attempted() == uninterrupted.attempted(),
        "final counters match the uninterrupted run");
}

void tamper_rejection(const std::string& path) {
  dqmc::Solver solver(test_params());
  for (int i = 0; i < 3; ++i) solver.sweep();
  solver.save_checkpoint(path, 3);

  std::ifstream in(path, std::ios::binary | std::ios::ate);
  const std::streamsize size = in.tellg();
  in.seekg(0);
  std::string bytes(static_cast<size_t>(size), '\0');
  in.read(bytes.data(), size);
  in.close();

  bytes[static_cast<size_t>(size) / 2] ^= 0x01;
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write(bytes.data(), size);
  }
  bool threw = false;
  try {
    dqmc::Solver::load_checkpoint(path);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "flipped payload byte is rejected");

  bytes[0] = 'X';
  {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write(bytes.data(), size);
  }
  threw = false;
  try {
    dqmc::Solver::load_checkpoint(path);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "bad magic is rejected");
}

void save_guard(const std::string& path) {
  dqmc::Parameters params = test_params();
  params.recompute_each_slice = false;
  params.recompute_after_sweep = false;
  dqmc::Solver solver(params);
  bool threw = false;
  try {
    solver.save_checkpoint(path, 0);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "save without recompute flags is refused");
}

}  // namespace

int main() {
  const std::string path = "test_checkpoint.bin";
  bitwise_trajectory(path);
  tamper_rejection(path);
  save_guard(path);
  std::remove(path.c_str());
  std::remove((path + ".tmp").c_str());
  if (failures == 0) {
    std::cout << "all checkpoint tests passed\n";
    return 0;
  }
  std::cout << failures << " checkpoint tests failed\n";
  return 1;
}
