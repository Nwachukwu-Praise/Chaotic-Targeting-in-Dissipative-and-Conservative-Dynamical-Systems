# Shinbrot Method Repair Audit

Date: 2026-08-14

This audit records the repair needed to separate a paper-faithful Shinbrot
bisection path from later project extensions. The primary authority used here
is `PhysRevA.45.4165.pdf`.

## Source Availability

- Reviewed: `PhysRevA.45.4165.pdf`.
- Reviewed current MATLAB files in `C:\Users\npcoz\OneDrive - Imperial College London\Project\Codes\Shinbrot targeting`.
- Not found locally: `Targeting__Overview_ (4).pdf`.
- Not found locally: `recovered_pre_august/Chaos report/run_shinbrot_discontinuity_aware_bisection(2).m`.
- Not found locally: `recovered_pre_august/Chaos report/search_parameter_to_target(8).m`.
- Not found locally: `recovered_pre_august/Chaos report/main_shinbrot_lorenz_demo(8).m`.

Statements below that refer to `Targeting__Overview_ (4).pdf` are therefore
requirements to check and correct in the report when that file is available,
not a direct page-by-page audit of an unavailable PDF.

## Paper-Grounded Corrections

- The paper applies the control parameter as an additive term in the Lorenz
  `Y` equation. It is not a state kick on a Poincare section.
- The paper section and target are `Z = 26.921`, accepted half-plane
  `X > 8.0`, target `X_t = 13.729`, and tolerance `0.008`.
- The published noise protocol retargets after every 40 integration time
  steps. It does not say 40 accepted Poincare crossings.
- The paper reports 24 refinements of the parameter range. The repaired
  paper path records the first tolerance hit separately and continues the
  final bracket refinement unless an effectively exact root is encountered.
- The paper detects discontinuities by comparing the ordered sequence of
  left and right half-plane crossings. Later whole-interval coverage,
  Lipschitz, dense probing, and hybrid-grid procedures must not be described
  as the original Shinbrot method.
- Raw MATLAB `sigmaNoise` is a per-step, per-coordinate coordinate standard
  deviation. It is not the normalized abscissa used in Shinbrot et al.
  Figure 6.

## File-By-File Audit

- `run_shinbrot_paper_bisection.m`: added as the separate paper-faithful
  deterministic path. It performs endpoint checks, crossing-order branch
  comparison, discontinuity range reduction, 24 final bracket refinements,
  and replay verification from the original source state.
- `search_parameter_to_target.m`: added dispatcher identifier
  `shinbrotPaperBisection`. The older `shinbrotDiscontinuityAwareBisection`
  remains available as a distinct project extension.
- `verified_bisection_identifier.m`: now returns `shinbrotPaperBisection`,
  so old saved "verified bisection" noise results are rejected unless they
  match the repaired implementation fingerprint and protocol.
- `run_noisy_targeting_trial.m`: added `paper40StepRetargetMode`. This mode
  uses a global RK4 step counter, adds Gaussian coordinate noise before each
  fixed RK4 step, detects section crossings by interpolation across the
  integrated step, and retargets on steps `40, 80, 120, ...`.
- `run_noise_sweep.m`: updated sweep validation, preallocated trial schema,
  trial tables, and provenance checks to support both
  `paper40StepRetargetMode` and the explicitly labelled
  `sectionRetargetMode` extension.
- `build_bisection_noise_config.m`: updated fingerprints to include mode
  list, RK4 step count, and 40-step retarget interval. This rejects stale
  checkpoints/results generated under the old section-retarget protocol.
- `load_valid_bisection_noise_results.m`: updated saved-result validation
  to require the current mode list and implementation fingerprint.
- `main_shinbrot_noise_demo.m`: default full and smoke sweeps now use
  `paper40StepRetargetMode` and new output filenames, preserving old
  section-retarget files for comparison.
- `plot_noise_results.m`: label-only update so selected-parameter plots do
  not imply section-only retargeting.
- `test_shinbrot_paper_repair.m`: added focused tests for dispatcher
  provenance, absence of forbidden paper-path machinery, replay verification,
  stale result rejection, reproducible independent Gaussian draws, and the
  40-step retarget schedule in a reduced smoke trial.

## Report Statements To Correct When The Report Is Available

- Replace any statement that the noise protocol retargets after "40 accepted
  Poincare crossings" with "40 integration time steps."
- Do not claim that the previous section-retargeting experiment reproduces
  the paper's noise-recovery experiment. It is a project extension.
- Do not describe whole-interval coverage, non-monotone feature subdivision,
  dense probes, Lipschitz bounds, or hybrid-grid searches as the original
  Shinbrot procedure.
- Do not compare raw `sigmaNoise` directly with the paper's normalized
  Figure 6 noise axis.
- Mark numerical tables and figures produced under old implementation
  fingerprints or the section-retarget protocol as stale for paper-replication
  claims.
- Do not describe replay through the same deterministic evaluator as fully
  independent validation. It is a reproducibility check under the same
  numerical model.
- Check and repair broken subsection and bibliography references.

## Remaining Numerical Limitations

- The deterministic search still uses `ode45` for high-accuracy verification.
  The paper used fixed-step RK4 and linear interpolation, but did not publish
  the RK4 step size.
- The crossing direction is retained as a project convention. The paper
  specifies the section and half-plane, but the provided text does not
  independently establish an upward-only crossing condition.
- Discontinuity range reduction needs a tie-breaking rule when both halves
  remain plausible. The repaired implementation documents its stack and
  closer-endpoint ordering as a project-specific choice, not as a published
  Shinbrot rule.
