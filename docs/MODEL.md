# Model and numerical conventions

## Target

Fontenele, Costa, dos Santos, and Paiva, Phys. Rev. B 105, 184502 (2022),
at `<n>=0.5`, `U/t=-5` on the periodic square lattice. Published values:
`Tc/t=0.150±0.003` from the helicity-modulus analysis in Fig. 3 and
`Tc/t=0.152±0.002` in Table 1.

## Hamiltonian

With `t=1` as the energy unit,

```
H = -t Σ_<ij>,σ (c†_iσ c_jσ + c†_jσ c_iσ)
    - μ Σ_iσ n_iσ
    - |U| Σ_i (n_i↑ - 1/2)(n_i↓ - 1/2).
```

The kinetic single-particle matrix is denoted `K`: off-diagonal entries are
`-t` on nearest-neighbor bonds and diagonal entries are `-μ`. Thus
`exp(-Δτ K)` contains the hopping and chemical-potential factors.

## Charge-channel discrete HS transformation

Let `a=Δτ |U|` and `cosh(λ)=exp(a/2)`. For every onsite occupation
`n↑,n↓ ∈ {0,1}`, the following exact identity holds:

```
exp[a (n↑-1/2)(n↓-1/2)]
  = (1/2) exp(-a/4) Σ_s=±1 exp[λ s (n↑+n↓-1)].
```

The four states check directly: identical occupations give `cosh(λ)` on the
right, while different occupations give unity. The factor `exp(-a/4)` drops
out of Metropolis ratios. The remaining `exp(-λ s)` normalization is
field-dependent and must be retained in the acceptance ratio.

For each imaginary-time slice `l`, define

```
B_l = exp(-Δτ K) D_l,
D_l(i,i) = exp(λ s_li).
```

The same matrix and determinant occur for each spin, so the total fermion
weight is the square of a real determinant. Define the cyclic product and
Green function at slice `l` by

```
C_l     = B_(l-1) … B_1 B_M … B_(l+1),
G_l     = (I + C_l B_l)^(-1).
```

For one noninteracting level this gives `G=1-f`, so onsite density is
measured as `n_iσ=1-G_l(i,i)`.

## Local Ising flip

At site `i` on slice `l`, let `s=s_li`, `r=exp(-2λ s)`, and `δ=r-1`.
Flipping `s -> -s` has the per-spin determinant ratio

```
R = 1 + δ [1-G_l(i,i)].
```

The determinant part of the two-spin Metropolis ratio is `R²`. Including the
HS normalization, the full ratio is `exp(2λs) R²`; this is the acceptance
weight used in the code. An accepted flip updates the Green function without
a full determinant recomputation:

```
G'_jk = G_jk - (δ/R) (δ_ji - G_ji) G_ik.
```

Advancing one slice uses

```
G_(l+1) = B_l G_l B_l^(-1),
```

with the same cyclic convention after wrapping at `l=M`. Periodic full
recomputation is a numerical-stability and correctness check, not an optional
substitute for testing the update formulas.

## Basic observables

- Density: `n=(1/N) Σ_iσ (1-G_ii)`.
- Double occupancy: onsite normal-ordered two-spin expectations must be
  evaluated with the interacting/tiny-exact reference, not inferred from the
  mean-field factorization.
- Pairing: `b_i=c_i↓ c_i↑`; pair correlations and `P_s` follow the paper's
  definitions in Eqs. (3)–(6). Signs at coincident sites are covered by tests.

## Superfluid stiffness

The target convention is

```
ρ_s = (1/4) [<-K_x> - Λ_T],
```

where `<-K_x>` is the x-bond kinetic (diamagnetic) term and `Λ_T` is the
zero-Matsubara-frequency transverse current-current response at `q_x=0`,
`q_y→0`. The independent longitudinal gauge identity at nonzero momentum is

```
Λ_L = <-K_x>.
```

The paramagnetic bond current is

```
j_x(l) = i t Σ_σ (c†_(l+x)σ c_lσ - c†_lσ c_(l+x)σ).
```

For momenta with `<J>=0`, full and connected current-current response coincide.
The implementation evaluates the connected Wick contraction

```
Λ = (2 Δτ/N) Σ_l -tr[J G(τ_l,0) J† G(0,τ_l)].
```

The time-displaced Green functions use

```
P_l = B_(l-1) … B_0
Q_l = B_(M-1) … B_l
G(τ_l,0) = (P_l^-1 + Q_l)^-1
G(0,τ_l) = -(Q_l^-1 + P_l)^-1.
```

`P` is accumulated with incremental LDR factorization. Each `Q` is independently
constructed in the required cyclic order with the same prepend operation. This
avoids forming or directly inverting a long, ill-conditioned product. Equal-time
Green functions use the corresponding cyclic LDR product and
`inverse_identity_plus`, which computes `(I+A)^{-1}`.

`Parameters::recompute_each_slice` is enabled by default. After each propagation
step, the Green function is recomputed from the full stabilized product. This
prevents update drift in long-slice runs; disabling it requires rerunning the
response and convergence tests.

The Kosterlitz–Thouless extraction locates the intersection of `ρ_s(T)` with
`2T/π`. A pair-structure-factor finite-size estimate is an independent check
when statistically usable.

## Validation status

The current tests cover the four-state HS identity, analytic free-fermion
checks, local determinant/Green updates against full recomputation on random
small fields, and tiny interacting-system observables against exact
diagonalization. Direct inverse and 40-factor LDR tests also pass, together with
exhaustive-field response through M=5, nonzero-longitudinal gauge, 2D transverse,
and independently nonzero 2D momentum checks. NaN and Inf are rejected in every
numerical comparison.

Checkpoint/restart is implemented: `src/checkpoint.cpp` writes versioned,
checksummed state, and `tests/test_checkpoint.cpp` verifies bitwise
trajectory equivalence between a resumed run and an uninterrupted one.
Production runs remain to be done at `<n>=0.5`, `U/t=-5`, `L=14,16,18`, and `Δτ=0.1`, with a documented
`Δτ=0.05` subset and independent seeds. Statistical, finite-size, and Trotter
sensitivity must accompany any comparison with the published `Tc` values.
