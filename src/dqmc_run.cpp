#include "dqmc/dqmc.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
  int lx = 0;
  int ly = 0;
  double temperature = -1.0;
  double beta = -1.0;
  double dtau = 0.1;
  double abs_u = 5.0;
  double mu = 0.0;
  bool mu_given = false;
  double density = -1.0;
  int warmup = 1000;
  int sweeps = 0;
  std::uint64_t seed = 12345;
  int measure_every = 1;
  int checkpoint_every = 500;
  std::string checkpoint_path;
  std::string resume_path;
  std::string out_path;
};

[[noreturn]] void usage(const char* argv0) {
  std::fprintf(stderr,
               "usage: %s --lx N --ly N (--T temp | --beta beta) --sweeps N\n"
               "  [--dtau 0.1] [--U 5] [--mu mu | --density n]\n"
               "  [--warmup 1000] [--seed 12345] [--measure-every 1]\n"
               "  [--checkpoint FILE [--checkpoint-every 500]] [--resume FILE]\n"
               "  [--out FILE]\n"
               "resume carries lattice, beta, dtau, U, mu, and seed from the "
               "checkpoint; --sweeps counts additional sweeps\n",
               argv0);
  std::exit(2);
}

Options parse(int argc, char** argv) {
  Options opt;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) usage(argv[0]);
      return argv[++i];
    };
    if (arg == "--lx") opt.lx = std::stoi(next());
    else if (arg == "--ly") opt.ly = std::stoi(next());
    else if (arg == "--T") opt.temperature = std::stod(next());
    else if (arg == "--beta") opt.beta = std::stod(next());
    else if (arg == "--dtau") opt.dtau = std::stod(next());
    else if (arg == "--U") opt.abs_u = std::stod(next());
    else if (arg == "--mu") { opt.mu = std::stod(next()); opt.mu_given = true; }
    else if (arg == "--density") opt.density = std::stod(next());
    else if (arg == "--warmup") opt.warmup = std::stoi(next());
    else if (arg == "--sweeps") opt.sweeps = std::stoi(next());
    else if (arg == "--seed") opt.seed = std::stoull(next());
    else if (arg == "--measure-every") opt.measure_every = std::stoi(next());
    else if (arg == "--checkpoint-every") opt.checkpoint_every = std::stoi(next());
    else if (arg == "--checkpoint") opt.checkpoint_path = next();
    else if (arg == "--resume") opt.resume_path = next();
    else if (arg == "--out") opt.out_path = next();
    else usage(argv[0]);
  }
  if (opt.sweeps <= 0 || opt.measure_every <= 0 || opt.checkpoint_every <= 0)
    usage(argv[0]);
  if (opt.resume_path.empty()) {
    if (opt.lx <= 0 || opt.ly <= 0) usage(argv[0]);
    if (opt.temperature <= 0.0 && opt.beta <= 0.0) usage(argv[0]);
    if (!opt.mu_given && opt.density <= 0.0) usage(argv[0]);
  }
  return opt;
}

struct Measurement {
  double density;
  double energy;
  double double_occupancy;
  double pair_structure;
  double minus_kinetic_x;
  double lambda_T;
  double lambda_L;
};

Measurement measure(const dqmc::Solver& solver) {
  Measurement m;
  m.density = solver.density();
  m.energy = solver.energy();
  m.double_occupancy = solver.double_occupancy();
  m.pair_structure = solver.pair_structure();
  m.minus_kinetic_x = solver.minus_kinetic_x();
  m.lambda_T = solver.current_response(0, 1).response;
  m.lambda_L = solver.current_response(1, 0).response;
  return m;
}

// rho_s = (<-K_x> - Lambda_T)/4 takes the TRANSVERSE response at qx=0,
// qy->0, so lambda_T_qy1 is the production channel. lambda_L_qx1 is the
// longitudinal response; the gauge identity Lambda_L = <-K_x> makes it a
// correctness check, not a production observable.
void write_record(std::ostream& out, std::uint64_t sweep,
                  const Measurement& m) {
  out << "{\"sweep\":" << sweep
      << ",\"density\":" << m.density
      << ",\"energy\":" << m.energy
      << ",\"double_occupancy\":" << m.double_occupancy
      << ",\"pair_structure\":" << m.pair_structure
      << ",\"minus_kinetic_x\":" << m.minus_kinetic_x
      << ",\"lambda_T_qy1\":" << m.lambda_T
      << ",\"lambda_L_qx1\":" << m.lambda_L << "}\n";
}

double stderr_of_mean(const std::vector<double>& xs, double mean) {
  if (xs.size() < 2) return 0.0;
  double var = 0.0;
  for (const double x : xs) var += (x - mean) * (x - mean);
  var /= static_cast<double>(xs.size() - 1);
  return std::sqrt(var / static_cast<double>(xs.size()));
}

struct Binning {
  double mean = 0.0;
  double stderr_naive = 0.0;
  double stderr_binned = 0.0;
  int bin_size = 1;
};

// Pairwise blocking with doubling bin size. The plateau is the standard
// error at the coarsest level still holding at least 8 bins, where the
// bin-to-bin autocorrelation is presumed small.
Binning bin_series(const std::vector<double>& series) {
  Binning result;
  if (series.empty()) return result;
  double sum = 0.0;
  for (const double x : series) sum += x;
  result.mean = sum / static_cast<double>(series.size());
  std::vector<double> level = series;
  result.stderr_naive = stderr_of_mean(level, result.mean);
  result.stderr_binned = result.stderr_naive;
  while (level.size() / 2 >= 8) {
    std::vector<double> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i + 1 < level.size(); i += 2)
      next.push_back(0.5 * (level[i] + level[i + 1]));
    level = std::move(next);
    result.bin_size *= 2;
    result.stderr_binned = stderr_of_mean(level, result.mean);
  }
  return result;
}

void write_summary_entry(std::ostream& out, const char* name,
                         const std::vector<double>& series) {
  const Binning b = bin_series(series);
  out << ",\"" << name << "\":{\"mean\":" << b.mean
      << ",\"stderr_naive\":" << b.stderr_naive
      << ",\"stderr_binned\":" << b.stderr_binned
      << ",\"bin_size\":" << b.bin_size << "}";
}

}  // namespace

int run(const Options& opt, std::chrono::steady_clock::time_point wall_start);

int main(int argc, char** argv) {
  const Options opt = parse(argc, argv);
  const auto wall_start = std::chrono::steady_clock::now();
  try {
    return run(opt, wall_start);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}

int run(const Options& opt, std::chrono::steady_clock::time_point wall_start) {
  std::optional<dqmc::Solver> solver;
  std::uint64_t start_sweep = 0;

  if (!opt.resume_path.empty()) {
    solver = dqmc::Solver::load_checkpoint(opt.resume_path);
    start_sweep = solver->sweep_tag();
  } else {
    const double beta =
        opt.beta > 0.0 ? opt.beta : 1.0 / opt.temperature;
    dqmc::Parameters params;
    params.lx = opt.lx;
    params.ly = opt.ly;
    params.dtau = opt.dtau;
    params.abs_u = opt.abs_u;
    params.seed = opt.seed;
    params.slices = static_cast<int>(beta / opt.dtau + 0.5);
    if (params.slices <= 0)
      throw std::invalid_argument("beta/dtau gives zero slices");

    if (!opt.mu_given) {
      const dqmc::CalibrationResult calibration = dqmc::calibrate_mu(
          params, opt.density, 100, 100, 1e-3, 20);
      params.mu = calibration.mu;
      std::fprintf(stderr, "calibrated mu = %.8f (density %.6f, %d iterations)\n",
                   calibration.mu, calibration.density, calibration.iterations);
    } else {
      params.mu = opt.mu;
    }
    solver.emplace(params);
  }

  const dqmc::Parameters& params = solver->params();
  const double beta = params.slices * params.dtau;

  std::ofstream out_file;
  std::ostream* out = &std::cout;
  if (!opt.out_path.empty()) {
    out_file.open(opt.out_path, std::ios::app);
    if (!out_file) throw std::runtime_error("cannot open output file");
    out = &out_file;
  }
  *out << "{\"header\":{\"lx\":" << params.lx << ",\"ly\":" << params.ly
       << ",\"slices\":" << params.slices << ",\"dtau\":" << params.dtau
       << ",\"beta\":" << beta << ",\"U\":" << params.abs_u
       << ",\"mu\":" << params.mu << ",\"seed\":" << params.seed
       << ",\"warmup\":" << opt.warmup
       << ",\"start_sweep\":" << start_sweep << "}}\n";

  for (int i = 0; i < opt.warmup && start_sweep == 0; ++i) solver->sweep();

  std::vector<double> density_series;
  std::vector<double> minus_kx_series;
  std::vector<double> lambda_T_series;
  std::vector<double> lambda_L_series;

  for (int i = 0; i < opt.sweeps; ++i) {
    solver->sweep();
    const std::uint64_t sweep = start_sweep + static_cast<std::uint64_t>(i) + 1;
    if (i % opt.measure_every == 0) {
      const Measurement m = measure(*solver);
      write_record(*out, sweep, m);
      out->flush();
      density_series.push_back(m.density);
      minus_kx_series.push_back(m.minus_kinetic_x);
      lambda_T_series.push_back(m.lambda_T);
      lambda_L_series.push_back(m.lambda_L);
    }
    if (!opt.checkpoint_path.empty() &&
        (i + 1) % opt.checkpoint_every == 0)
      solver->save_checkpoint(opt.checkpoint_path, sweep);
  }
  if (!opt.checkpoint_path.empty())
    solver->save_checkpoint(opt.checkpoint_path,
                            start_sweep + static_cast<std::uint64_t>(opt.sweeps));

  if (!density_series.empty()) {
    *out << "{\"summary\":{\"measurements\":" << density_series.size();
    write_summary_entry(*out, "density", density_series);
    write_summary_entry(*out, "minus_kinetic_x", minus_kx_series);
    write_summary_entry(*out, "lambda_T_qy1", lambda_T_series);
    write_summary_entry(*out, "lambda_L_qx1", lambda_L_series);
    *out << "}}\n";
    out->flush();
  }

  const double wall = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - wall_start).count();
  std::fprintf(stderr,
               "done: %llu sweeps in %.1f s (%.3f s/sweep), acceptance %.4f\n",
               static_cast<unsigned long long>(opt.sweeps), wall,
               wall / opt.sweeps,
               solver->attempted() > 0
                   ? static_cast<double>(solver->accepted()) /
                         solver->attempted()
                   : 0.0);
  return 0;
}
