# ECR-III for the Lorenz system — MATLAB implementation

A from-scratch MATLAB (and GNU Octave) implementation of the **ECR-III**
targeting method of

> S. Iplikci, Y. Denizhan, *An improved neural network based targeting method
> for chaotic dynamics*, Chaos, Solitons and Fractals **17** (2003) 523–529,

applied to the **Lorenz system**, in real coordinates and in delay coordinates.
No toolboxes are required: the RBF networks, the k-means, the clustering and the
Poincaré map are all in this folder.

ECR-III in one paragraph: the phase space is carved into *control regions*
`T_0, T_1, …, T_K`, where `T_0` is the small target region around the unstable
orbit and `T_i` is the set of states that can be pushed into `T_{i-1}` in one
step with an admissible parameter change. Each region is split into *clusters*
by a simple radius rule, each cluster gets its own small RBF network mapping
state → parameter change, and each cluster is summarised by its mean and its
normalised covariance. On line, the current state is assigned to the cluster with
the smallest **Normalised Mahalanobis Distance** (Eq. 8) and that cluster's
network supplies the parameter perturbation. Nothing but recorded input/output
data is used — no model, no linearisation.

---

## Quick start

```matlab
>> cd ECR_III_MATLAB
>> startup_ecr                 % puts core/ systems/ demos/ tests/ on the path
>> demo_lorenz_ecr3            % the main experiment (~5 min)
```

Other entry points:

| script | what it does |
|---|---|
| `ECR_III_walkthrough` | narrative walkthrough of the whole method with equations and inline figures — the live-script source (see below) |
| `demo_lorenz_ecr3` | full ECR-III experiment in real coordinates: data → regions → clusters → networks → reaching times → figures |
| `demo_lorenz_delay_ecr3` | the same with a single measured signal and a 100 ms delay embedding |
| `run_lorenz_comparison` | Table-1 style sweep over the number of control regions, ECR-II vs ECR-III, noiseless and noisy → `results/*.csv` |
| `demo_options_study` | measures the effect of the choices the paper leaves open (see `docs/METHOD_NOTES.md`) |
| `test_all` | the test suite (`test_all(true)` for the quick version) |

Generated data are cached in `results/data_lorenz_real.mat`; delete the file to
regenerate.

### Live script

`demos/ECR_III_walkthrough.m` is written in live-script markup (section
headings, formatted text, LaTeX equations, inline figures). A `.mlx` is a
binary package that only MATLAB itself can write, so convert it once:

```matlab
>> startup_ecr
>> make_livescript                 % demos/ECR_III_walkthrough.m -> .mlx
>> open ECR_III_walkthrough.mlx
```

`make_livescript` calls the Live Editor conversion API; if your release does
not expose it, open the `.m` in MATLAB and use **Save As → MATLAB Live Code
File (\*.mlx)** — same result. The `.m` also runs as an ordinary script in
MATLAB and Octave.

---

## What you get

```
ECR_III_MATLAB/
├── startup_ecr.m              add everything to the path
├── systems/
│   ├── sys_lorenz_real.m      Lorenz as a Poincare map on y = 8.4853, z = [x;z]
│   ├── sys_lorenz_delay.m     Lorenz seen through one scalar signal, delay embedded
│   ├── lorenz_ode.m / _vec.m  the flow (scalar and vectorised)
│   ├── flow_rk4_step*.m       fixed-step RK4
│   ├── flow_poincare_step*.m  integrate to the next section crossing (Newton-refined)
│   └── ecr_system_template.m  the struct every system fills in
├── core/
│   ├── ecr_default_options.m  every tunable knob, with the interpretation flags
│   ├── ecr_generate_data.m    random-parameter experiments -> (z_n, p_n, z_{n+1})
│   ├── ecr_train.m            Definitions 1 and 3: builds T_0 … T_K and their nets
│   ├── ecr_choose_radius.m    clustering radius from the inter-distance histogram
│   ├── ecr_cluster.m          the paper's simple r-neighbourhood clustering
│   ├── ecr_cluster_stats.m    cluster means and normalised covariance matrices
│   ├── ecr_nmd.m              Eq. (8)
│   ├── ecr_level_eval.m       "does this region recognise the state?" (Fig. 3)
│   ├── ecr_control.m          the on-line controller (Fig. 4), ECR-II and ECR-III
│   ├── ecr_simulate.m         one closed-loop run
│   ├── ecr_reaching_time.m    average reaching time benchmark (Section 4)
│   ├── rbf_train.m / rbf_eval.m / rbf_design.m / kmeans_simple.m
│   ├── ecr_plot_clusters.m / ecr_plot_run.m / ecr_summary.m
│   └── ecr_cache_data.m / ecr_save_model.m / ecr_load_model.m / ecr_seed.m
├── demos/                     the walkthrough + the four scripts listed above
├── tools/make_livescript.m    convert the walkthrough into a .mlx live script
├── tests/                     test_all.m and seven test files (+ the logistic fixture)
├── docs/METHOD_NOTES.md       equation-by-equation mapping + every judgement call
└── results/                   figures, csv tables and cached data (created on the fly)
```

---

## Reproduction status

| quantity | paper (Table 1 / 2) | this code |
|---|---|---|
| Poincaré-map target `z*` | `[14.2387 39.7934]` | `[14.2522 39.7867]` |
| RMS of the measured state | 13.43 | ≈13.5 |
| reaching time, no targeting | 15.09 (OGY column) | ≈15–20 |
| reaching time, ECR-II | 5.9–6.7 | 4.5–6.0 |
| reaching time, ECR-III | 6.4–6.8 | 4.1–5.3 |
| delay coordinates: none / ECR-II / ECR-III | 35.3 / 12.8–18.9 / 12.5–16.2 | ≈11 / 7.8 / 5.6 |
| logistic-map fixture, ECR-III / ECR-II / none | 5.52 / 9.14 / 143.7 | ≈5.7 / ≈10.3 / ≈105 |

The delay-coordinate target of Table 2 could **not** be reproduced: the paper
does not say which signal was measured, and no coordinate of the Lorenz system
matches the quoted values on the stated surface of section. `sys_lorenz_delay`
keeps the structure and computes its own target — details in
`docs/METHOD_NOTES.md`, which also lists every place where the paper is silent
and what this implementation does instead.

## Notes

* Out of scope by request: OGY control, ECR-I, and every system other than the
  Lorenz one. Inside the target region the local controller is `NN_0`, trained
  from data exactly like the other region networks. The logistic map remains
  only as a fast test fixture.
* Everything runs in MATLAB and in Octave ≥ 7 (the suite was developed and
  verified against Octave 8.4).
* A full `demo_lorenz_ecr3` run is a few minutes; the delay-coordinate demo is
  slower because every map step integrates the flow while recording a
  measurement history.
* In delay coordinates ECR-III reaches the target quickly but does not *hold*
  it: a delay vector is partly generated under the previously applied
  parameter, so the observed map is not a function of the current state and
  parameter alone. This is a property of delay-coordinate parameter control,
  not of ECR — see the last section of `docs/METHOD_NOTES.md`.
