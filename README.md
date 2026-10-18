# attractive-hubbard-dqmc

Serial determinant quantum Monte Carlo for the two-dimensional attractive
Hubbard model. The initial target is Fontenele et al., Phys. Rev. B 105,
184502 (2022), at density `<n>=0.5` and `U/t=-5`, where Fig. 3 reports
`Tc/t=0.150±0.003` and Table 1 reports `Tc/t=0.152±0.002`.

This code is independent of the harmonic-oscillator DMC payload in `qmc-rdma`.
GPU and RDMA acceleration are not part of the current serial implementation.

## Current status

The determinant engine, exact references, pairing observables,
chemical-potential calibration, and current-response implementation have passed
the existing numerical checks. CTest passes 2/2.

- Shared incremental LDR factorization and direct inverse tests cover identity,
  diagonal, permutation, well-conditioned, and 40-factor product cases.
- The response harness rejects NaN/Inf, checks nonzero longitudinal momentum,
  and includes a 2x2 transverse case plus an independently nonzero 2D momentum
  response.
- Exhaustive-field results through M=5 and long-slice scans are finite. At
  beta=8, the response converges `1.38021 -> 1.29356 -> 1.28083` toward the
  exact value `1.27613` for dtau `0.2 -> 0.1 -> 0.05`.

The code does not yet provide a production helicity-modulus or `Tc` result,
checkpoint/restart equivalence, or a comparison with the published `Tc` values.
See `HANDOFF.md` for the validation record and remaining production work, and
`docs/MODEL.md` for the equations and conventions.

## Build and test

The implementation uses C++20 and Eigen, with CMake/CTest and Python
standard-library-only analysis scripts.

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

Green recomputation after every slice is enabled by default to prevent update
drift in long-slice runs. It is slower than update-only propagation and should
not be disabled without rerunning the response and convergence tests.
