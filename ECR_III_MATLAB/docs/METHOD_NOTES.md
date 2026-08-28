# ECR-III implementation notes

Reference: S. Iplikci, Y. Denizhan, *"An improved neural network based targeting
method for chaotic dynamics"*, Chaos, Solitons and Fractals **17** (2003) 523–529.

This document maps every element of the paper onto the MATLAB code, and lists —
explicitly — every place where the paper is silent and a decision had to be made.
Nothing here is hidden inside the code: each decision is an option in
`ecr_default_options.m` and is marked `[interpretation]` there.

---

## 1. What is implemented

| Paper | Code |
|---|---|
| Eq. (1) `z_{n+1} = G(z_n,p_n)` | `S.step(H,P)` in `systems/sys_lorenz_*.m` |
| Eq. (2) target `z*` | `S.zstar`, found by Newton on the Poincaré map |
| Eq. (3) allowable set `Π` | `ecr_in_pi.m`, `S.pnom`, `S.dpmax` |
| Eq. (4) Definition 1, region `S0 = T0` | level-0 selection in `ecr_train.m` |
| Eq. (7) Definition 3, regions `T_i` | level loop in `ecr_train.m` |
| Fig. 3 data extraction with `NN_{i-1}` | `ecr_level_eval.m` used on the successor states |
| Section 3.3 radius from the inter-data-distance histogram | `ecr_choose_radius.m` |
| Section 3.3 "simple clustering algorithm" | `ecr_cluster.m` (r-neighbourhood chaining) |
| Section 3.3 cluster means / normalised variances | `ecr_cluster_stats.m` |
| Eq. (8) Normalised Mahalanobis Distance | `ecr_nmd.m` |
| Section 3 RBF networks `NN_ij` | `rbf_train.m`, `rbf_eval.m` (no toolbox) |
| Fig. 4 on-line control loop | `ecr_control.m`, `ecr_simulate.m` |
| Section 4 average reaching time | `ecr_reaching_time.m` |
| ECR-II (Section 3.2), used as the comparison method | `ecr_train(...,'ECR-II')` + `ecr_control` |

**Not implemented (out of scope by request):** the OGY controller and the ECR-I
single-network variant, and every system other than the Lorenz system. In the
paper the local controller inside `S0` is OGY; here `NN_0` — trained exactly like
every other region network, from data — plays that role, so no linearisation and
no model knowledge is used anywhere in the toolbox. The logistic map survives
only as a fast test fixture (`tests/fixtures/sys_logistic_fixture.m`), because a
1-D map makes the region logic easy to check and the paper quotes numbers for it.

---

## 2. The Lorenz set-up (Table 2 of the paper)

* Flow: `dx/dt = σ(y−x)`, `dy/dt = ρx − y − xz`, `dz/dt = xy − βz`.
* Control parameters `p = [σ ρ β]`, `p_nom = [10 28 8/3]`,
  `δp_max = [0.30 0.84 0.08]`, `δ = 0.30`.
* Surface of section `y = y_PSS = 8.4853 = sqrt(β(ρ−1))`.

**Crossing direction.** The paper does not state it. Both were tried:

* `sec_dir = −1` (crossings with `dy/dt < 0`) gives a Poincaré map whose fixed
  point is `z* = [14.2522 39.7867]`, versus `[14.2387 39.7934]` in Table 2, and
  an x-RMS of ≈13.5 versus 13.43 in Table 2. This is the default.
* `sec_dir = +1` puts the `C+` equilibrium `(8.4853, 8.4853, 27)` on the section
  itself, and the fixed-point search collapses onto it. It reproduces nothing in
  Table 2.

The target is therefore a period-1 unstable periodic orbit of the flow, not an
equilibrium (Table 2 calls the column "equilibrium point", but `[14.2387
39.7934]` is not an equilibrium of the Lorenz system). Its Poincaré-map
eigenvalues here are `λ ≈ 4.71` and `λ ≈ 10⁻⁹`: strongly unstable in one
direction, strongly contracting in the other.

**Integration.** Fixed-step RK4, `dt = 2·10⁻³` (real coordinates) or `10⁻³`
(delay coordinates), with the crossing time refined by Newton iterations on the
section function. Bisection was tried first and is *not* good enough: the
residual crossing-time error swamps finite-difference Jacobians of the map, which
made the fixed point look like a repeller in both directions. All trajectories of
a batch are integrated in lock-step (`flow_poincare_step_vec.m`), which is what
makes the experiments finish in minutes rather than hours.

**Delay coordinates.** Table 2 quotes `T = 100 ms`, an embedding
`z* = [−3.0988 −3.2911 −3.8340]` and a measured-signal RMS of 2.675. Those
numbers cannot be reproduced with `x`, `y` or `z` measured on the `y = 8.4853`
section: the delay vector of the period-1 orbit for `s = x` is
`[14.2522 13.3172 6.9024]`, and no coordinate of the Lorenz system has an RMS of
2.675 on that section. The paper does not say which signal was measured or how it
was scaled. `sys_lorenz_delay.m` therefore keeps the structure (scalar
measurement, `T = 100 ms`, embedding dimension 3) and computes its own target
from the same period-1 orbit; the measured signal is an option (`opts.meas`,
default `x`).

---

## 3. How the regions are built from data only

Data: `ecr_generate_data.m` runs the system from many initial conditions while a
parameter vector drawn uniformly from `Π` is applied at every step, and records
the triples `(z_n, p_n, z_{n+1})`. This is the only information about the plant
that the method ever uses.

* **T0 (Definition 1).** A triple qualifies if `‖z_n − z*‖ < δ` *and*
  `‖z_{n+1} − z*‖ < δ`. The recorded `p_n` is both the evidence that an
  admissible parameter exists and the training target for `NN_0`.
* **T_i (Definition 3).** A triple qualifies if the successor `z_{n+1}` is
  recognised by the already-trained level `i−1` model (its network returns a
  parameter inside `Π`) and `z_n` is not recognised by any lower level. This is
  exactly the data-extraction procedure of Fig. 3; no system equations are used.
* Every region is then clustered (ECR-III) and one RBF network is trained per
  cluster on `z ↦ (p − p_nom)/δp_max`, so "inside `Π`" is simply
  `max|output| ≤ 1`.

### Enrichment of the target region `[interpretation]`

A chaotic trajectory spends only a few percent of its time inside the `δ`-ball,
and only a fraction of *those* visits happen to receive a parameter that keeps it
there, so `T0` would be starved (≈30 triples out of 6000). `opt.enrich`
therefore re-runs states that came close to the target with fresh random
parameters, and — when the observed state can be lifted back to a full state
(real coordinates) — seeds extra experiments inside and around the `δ`-ball.
This is a *data-collection* convenience, not model knowledge: it only says "run
the experiment again from a similar state".

---

## 4. Decisions the paper leaves open

All of these are options; the defaults are the ones measured to work best on the
Lorenz system in `demo_options_study.m` (`results/lorenz_real_options.csv`).

| Option | Question the paper does not answer | Default | Why |
|---|---|---|---|
| `covNorm` | what "Normalized Covariance Matrix" means in Eq. (8) | `'det'` | `Σ/det(Σ)^{1/N}` removes the *size* of a cluster but keeps its *shape*, so clusters of different populations compete fairly; for an isotropic cluster the NMD reduces to the plain Euclidean distance. `'plain'` (ordinary Mahalanobis), `'diag'` ("normalised variances", no cross terms) and `'trace'` are available. |
| `nmdGate` | whether a cluster may be rejected because the state is far outside it | `2.0` | Eq. (8) always returns a nearest cluster, even for a state nowhere near any data; without a gate an RBF network can happily extrapolate a plausible-looking parameter into `Π`. The gate accepts a cluster only if the NMD is at most 2× the 95th percentile NMD of its own training data. `Inf` reproduces the literal paper. |
| `nFallback` | what to do when the nearest cluster's network returns `p ∉ Π` | `3` | The paper stops there ("no targeting is applied"). Trying the next-nearest clusters costs nothing and reduces the number of wasted steps. `1` = literal paper. |
| `selection` | which cluster wins when clusters of several regions are close | `'nmd'` | The literal Section 3.3 rule: minimum NMD over *all* clusters of *all* regions. `'level'` prefers the lowest region index (ECR-II ordering) and is a documented alternative. |
| `targetTighten` | which `S0` pairs train the local controller | `0.5` | Definition 1 only requires the successor to stay inside the `δ`-ball; pairs that land well inside it teach a much sharper local controller (retention after capture rises from ≈0.8 to ≈0.97). `1` = literal Definition 1. |
| `escapeAfter` | nothing — this is a benchmarking safeguard | `12` | Deterministic feedback can lock a trajectory into a *controlled* periodic orbit that never enters `S0`; the paper's tables show no such failures. After 12 consecutive controlled steps without capture the nominal parameters are applied for one step. `0` disables it. |
| clustering post-processing | what happens to clusters with a handful of points | merge into nearest | The paper's radius rule alone can leave dozens of one-point clusters, each producing a degenerate covariance and a constant network. Clusters below `minClusterPts` are merged into the nearest larger one; if *no* cluster is large enough the largest one absorbs the rest. |
| radius `r` | "first minimum of the histogram ... by visual analysis" | automated | The histogram of pairwise distances is smoothed with a 3-bin moving average and the first local minimum after the first local maximum is taken; `ecr_choose_radius` returns the histogram so the choice can still be inspected. If the histogram is single-moded the fallback is twice the median nearest-neighbour distance. |
| RBF details | centres, widths, weights | k-means + nearest-centre widths + ridge least squares | The paper only says "RBF-based neural networks ... better local approximation capability". Two-stage training (unsupervised centres, closed-form weights) is the standard choice and is what makes ECR training cheap, which is the property Table 1 measures. |

---

### Measured effect of those choices

`demo_options_study.m`, Lorenz in real coordinates, 6000 raw triples, five
control regions, 100 fixed initial conditions, each configuration trained three
times (average reaching time ± spread over the three trainings; "hold" is the
fraction of steps that stay inside the target region after capture):

| covNorm | gate 1.5 | gate 2.0 | gate ∞ | hold |
|---|---|---|---|---|
| `det`   | 4.57 ± 0.64 | **4.39 ± 0.43** | 4.77 ± 0.90 | 0.96–0.97 |
| `trace` | 4.46 ± 0.50 | **4.26 ± 0.23** | 4.50 ± 0.67 | 0.96 |
| `plain` | 5.41 ± 2.41 | 5.18 ± 1.30 | 5.26 ± 2.06 | 0.61–0.69 |
| `diag`  | 5.30 ± 1.94 | 5.17 ± 1.02 | 5.27 ± 1.87 | 0.67–0.68 |

(no targeting at all: 18.4 steps)

The two *size-removing* normalisations (`det`, `trace`) are both better and far
more repeatable than the two that keep the cluster size (`plain`, `diag`), and
they hold the orbit once captured. `det` and `trace` are statistically
indistinguishable here; `det` is kept as the default because a unit-determinant
covariance is the more natural reading of "normalised covariance matrix".

Region selection at the best setting: `'nmd'` — the literal Section 3.3 rule —
gives 4.26 ± 0.23 versus 4.81 ± 0.26 for the ECR-II-style `'level'` ordering.
The paper's rule wins, so it is the default.

## 5. What was reproduced

Using `demo_lorenz_ecr3.m` (6000 raw triples, 100 initial conditions,
`δ = 0.30`, five control regions):

* **Target.** `z* = [14.2522 39.7867]` versus Table 2's `[14.2387 39.7934]`;
  x-RMS on the section 13.5 versus 13.43.
* **Uncontrolled reaching time.** ≈15–20 steps depending on the sample of
  initial conditions, versus 15.09 for the OGY column of Table 1 (which is the
  same experiment: no targeting, wait until the trajectory enters the region).
* **ECR-II / ECR-III reaching time.** 4–7 steps, versus 5.9–6.8 in Table 1.
  As in the paper, for the Lorenz system in real coordinates the two methods are
  close — the Lorenz Poincaré map returns to the neighbourhood of the target
  quickly, so there is little room left between one region and six. The clear
  ECR-III advantage the paper reports for the logistic and Hénon maps is
  reproduced on the logistic fixture (see `tests/test_pipeline_logistic.m`:
  ≈5.7 steps for ECR-III versus ≈10.3 for ECR-II and ≈105 uncontrolled;
  the paper quotes 5.52 / 9.14 / 143.7).
* **Training cost.** ECR-III trains many small networks rather than one large
  one; with the closed-form RBF solution the whole model takes a few seconds.
  `results/lorenz_real_comparison.csv` (produced by `run_lorenz_comparison`)
  holds the sweep over the number of control regions:

  | regions | ECR-II reach | ECR-III reach | ECR-II train [s] | ECR-III train [s] |
  |---|---|---|---|---|
  | 1 | 5.95 | 5.31 | 0.11 | 0.46 |
  | 2 | 5.13 | 4.10 | 0.22 | 2.81 |
  | 3 | 4.54 | 5.14 | 0.23 | 3.75 |
  | 4 | 5.30 | 6.94 | 0.40 | 3.25 |
  | 5 | 5.93 | 5.30 | 0.41 | 4.23 |

  with 15.5 steps for no targeting. As in Table 1 the gain saturates after one
  or two regions for this system, and the two methods stay within about a step
  of each other — the Lorenz Poincaré map returns close to the target quickly,
  so there is little room for a deep region hierarchy. (In this implementation
  ECR-III does *not* train faster than ECR-II, unlike Table 1: the paper's
  ECR-II networks are trained iteratively, where cost grows steeply with the
  size of a region, while here every network — large or small — is fitted by
  the same closed-form least squares, so splitting a region into clusters
  cannot save time.)

* **Delay coordinates.** `demo_lorenz_delay_ecr3` gives ≈5.6 steps for ECR-III
  and ≈7.8 for ECR-II against ≈11 with no targeting, i.e. the ECR-III advantage
  the paper reports for this case is reproduced. The *retention* after capture,
  however, is poor (≈0.2 against ≈0.96 in real coordinates), and more data does
  not fix it. The reason is structural rather than a bug: a delay vector
  measured at step `n` was partly generated while the *previous* parameter
  value was applied, so the observed map `z_n -> z_{n+1}` is not a function of
  `z_n` and `p_n` alone. Coarse targeting survives that ambiguity, precise
  local control does not. The standard remedy is the Dressler–Nitsche
  correction (augment the observed state with the previously applied
  parameter); it is not implemented here because the paper does not use it.

## 6. Known limitations

* The reaching times depend on the random initial-condition sample and on the
  training seed; differences below roughly one step are not significant with 100
  trials. Fix `opt.seed` (and the seed used for `H0`) when comparing variants.
* `sys_lorenz_delay` reproduces the *structure* of the delay-coordinate
  experiment but not the exact numbers of Table 2 (see §2).
* Training-time comparisons are wall-clock, not the flop counts of Table 1.
* The clustering is the "simple" one of the paper; the authors themselves note
  that a better clustering algorithm would improve the reaching time further at
  the price of longer training.
