# Validation and reproduction notes

Updated 2026-10-27.

The determinant core and connected current-response path pass CTest 3/3.
Long-slice responses are finite, and the beta=8 response approaches the exact
reference as dtau decreases. Production helicity-modulus and `Tc` calculations
remain to be done.

## Reproduction target

Fontenele, Costa, dos Santos, and Paiva, Phys. Rev. B 105, 184502 (2022), report
results for the 2D square-lattice attractive Hubbard model at `<n>=0.5` and
`U/t=-5`:

- Fig. 3: `Tc/t=0.150±0.003`
- Table 1: `Tc/t=0.152±0.002`

The first production milestone is one critical-temperature point from the
helicity-modulus/Kosterlitz–Thouless analysis, not the full phase diagram.
GPU/RDMA acceleration and d-wave/cuprate interpretations are outside this target.

## Source layout

- `docs/MODEL.md`: model, HS transformation, Green/update conventions, and
  response equations
- `include/dqmc/dqmc.hpp`, `src/dqmc.cpp`: lattice, exact references,
  determinant engine, local updates, basic observables, pairing, calibration,
  and per-slice stabilized Green recomputation
- `include/dqmc/ldr.hpp`: shared one-shot and incremental LDR factorization and
  inverse helpers
- `include/dqmc/current.hpp`, `src/current.cpp`: connected current-current
  response, exact tiny Kubo reference, and deterministic test support
- `tests/test_core.cpp`: core correctness tests
- `tests/test_current.cpp`: direct inverse tests, exhaustive-field response,
  nonzero longitudinal momentum, and 2D response tests

## Core numerical validation

The core tests cover the charge-channel discrete HS identity, determinant
ratios, Green updates, free-fermion and tiny interacting references, chemical
potential, pair operators, pair structure, and chemical-potential bisection.

The current implementation also retains checks for the `exp(2λs) R²` acceptance
ratio, column scaling for `B=A D`, cyclic `A/A^{-1}` propagation order,
one-dimensional self-hopping, diagonal chemical-potential terms, two-site-ring
current tests, and inverse column permutation in one-shot LDR reconstruction.

## Inverse identities and stabilized products

- `inverse_identity_plus` computes `(I+A)^{-1}`. Deterministic tests compare it
  with direct inversion on identity, diagonal, permutation, and
  well-conditioned matrices. An earlier probe used `(A+A^2)^{-1}` and did not
  test the intended identity.
- `inverse_inverse_plus(u,v)` computes `(u^{-1}+v)^{-1}` and is tested on
  identity, diagonal, permuted, and well-conditioned pairs. A near-singular pair
  is not used as an absolute-error acceptance case.
- Incremental `prepend_ldr` folds each new slice matrix into the decomposition
  with a fresh QR instead of forming an ill-conditioned long product. A
  40-factor reconstruction test passes at `1e-10` tolerance.
- `green_at_slice` uses the same stabilized decomposition. Green recomputation
  after every slice is enabled by default through
  `Parameters::recompute_each_slice`; the previous update-only path drifted in
  long-slice runs.

## Connected current response

For nonzero momenta with `<J>=0`, the full and connected current-current
responses coincide. The implementation evaluates

```
Λ = (2 Δτ/N) Σ_l -tr[J G(τ_l,0) J† G(0,τ_l)].
```

The slice products are

```
P_l = B_(l-1) … B_0
Q_l = B_(M-1) … B_l
G(τ_l,0) = (P_l^-1 + Q_l)^-1
G(0,τ_l) = -(Q_l^-1 + P_l)^-1.
```

`P` is built incrementally. Each `Q` is built independently with the same
prepend operation to preserve cyclic order. No direct long-product inverse is
used. The exact Kubo reference is compared only where the exact current
expectation vanishes, so full and connected response are consistent.

The response tests:

- reject NaN and Inf in every numerical comparison;
- use `qx=1` for the longitudinal gauge identity rather than `q=0`;
- test a 2x2 transverse response and an independently nonzero 2D momentum
  response;
- check the 2D diamagnetic term against exact diagonalization.

Current CTest results:

- Longitudinal: exact `1.27068`, DQMC `1.29743`, kinetic `1.27704`.
- Nonzero 2D momentum: exact `0.987279`, DQMC `1.01129`, kinetic `0.986929`.
- The 2x2 transverse response is exactly zero by symmetry; DQMC also returns
  zero. The independently nonzero 2D check ensures that the test does not pass
  by returning zero everywhere.

## Enumeration and convergence

Exhaustive auxiliary-field enumeration gives:

| Slices | Enumerated response | Exact response |
|--:|--:|--:|
| 1 | 0.205574 | 0.197193 |
| 2 | 0.397493 | 0.379462 |
| 3 | 0.566857 | 0.538892 |
| 4 | 0.712462 | 0.674879 |
| 5 | 0.836922 | 0.790325 |

The independent sampling scan at dtau `0.2` gives finite responses through
M=40: M=3 `0.566272`, M=5 `0.843724`, M=10 `1.22185`, M=20 `1.35449`, and
M=40 `1.37841`. At fixed beta=8, reducing dtau gives:

| Slices | dtau | DQMC response | Exact response |
|--:|--:|--:|--:|
| 40 | 0.20 | 1.38021 | 1.27613 |
| 80 | 0.10 | 1.29356 | 1.27613 |
| 160 | 0.05 | 1.28083 | 1.27613 |

These are Trotter-convergence checks, not a production `Tc` result or a claim of
exact agreement at every Trotter step.

## Remaining production work

1. ~~Add and verify checkpoint/restart equivalence and a production CLI/data
   path.~~ **Done 2026-10-27.** `src/checkpoint.cpp` adds versioned,
   checksummed binary checkpoints (FNV-1a over the payload, atomic
   tmp+rename write). The stored state is parameters, Ising fields, the
   accumulated B matrices, the mt19937_64 stream, and the acceptance
   counters plus a driver sweep tag — the B matrices are stored because
   their incremental column updates accumulate rounding that a rebuild from
   fields would not reproduce, and bitwise continuation is the design
   requirement. `tests/test_checkpoint.cpp` verifies that a split run
   reproduces the uninterrupted observable stream exactly (density, energy,
   double occupancy, pair structure, -<K_x> after every sweep), that
   counters and the sweep tag round-trip, and that a flipped payload byte
   or bad magic is rejected. Checkpoints are valid at sweep boundaries and
   saving requires at least one recompute flag (the boundary Green function
   must be a pure function of the stored state). The production driver is
   `src/dqmc_run.cpp` (`build/dqmc_run`): lattice/beta-or-T/dtau/U/seed,
   `--mu` or `--density` (bisection calibration, coarse-grained on small
   lattices), warmup/sweeps/measure-every, `--checkpoint`/
   `--checkpoint-every`, `--resume` (sweeps count additionally), and JSONL
   records carrying the helicity-modulus ingredients (density, energy,
   double occupancy, pair structure, -<K_x>, and both connected current
   responses at the smallest nonzero momentum per measurement: the
   transverse Lambda_T at qx=0, qy=1 — the production channel entering
   rho_s = (<-K_x> - Lambda_T)/4 — and the longitudinal Lambda_L at qx=1,
   recorded as the gauge-identity check), closing with a blocking/binning
   summary of the measured series. Full suite: CTest 3/3.
2. ~~Measure performance without changing the validated numerical path.~~
   **Done 2026-10-27 (quiet machine, seconds/sweep, dtau=0.1, T=0.15,
   M=67, warmup excluded):** L=14 (N=196): 15.0; L=16 (N=256): 77.7;
   L=18 (N=324): 137.3. The dominant cost is the per-slice stabilized Green
   rebuild (correctness choice — do not disable without rerunning the
   response suite). L=14→16 is worse than N^3 (5.2x for 1.31x sites); the
   16→18 step is 1.8x — small-sample and thermal noise on a laptop, so treat
   the exponent as unresolved. Production cost at 10k sweeps per (T, seed)
   point: L=14 ~42 h, L=16 ~9 days, L=18 ~16 days CPU time. The temperature
   grid and seeds are independent single-core jobs, so wall time divides by
   core count; even so L=18 is not feasible at these sweep counts without a
   validated recompute-cadence reduction or a smaller grid.
3. Run `L=14,16,18` at `<n>=0.5`, `U/t=-5`, `Δτ=0.1`, with a documented
   `Δτ=0.05` subset and independent seeds. Grid and statistics sized from
   the phase-2 baseline above and an L=14 autocorrelation shakedown before
   any long launch.
4. Extract `ρ_s(T)` intersections with `2T/π`, reporting statistical,
   finite-size, and Trotter sensitivity. Compare pair-structure finite-size
   estimates when statistically usable.
5. Record the result, including negative or unresolved outcomes, and compare it
   with the published values without tuning temperatures, sweeps, or tolerances
   toward agreement.

## Build and test

```bash
cmake -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXE_LINKER_FLAGS=-fuse-ld=mold
cmake --build build
ctest --test-dir build --output-on-failure
```

The default per-slice recomputation is slower than update-only propagation. Do
not disable it without rerunning the response and convergence tests.
