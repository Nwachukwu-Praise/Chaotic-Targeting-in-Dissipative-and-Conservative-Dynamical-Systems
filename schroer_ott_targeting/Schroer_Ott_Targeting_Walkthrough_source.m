%% Schroer-Ott targeting for the lifted standard map
% This Live Script rebuilds the corrected Schroer-Ott workflow from the
% computational pieces rather than calling |main_schroer_ott_targeting|.
%
% The state space used here is the lifted cylinder S^1 x R: x is periodic
% modulo one, while y is never wrapped.  Controls are instantaneous vertical
% y-kicks with |u| <= 0.003.  The operational targeting method uses small
% proxy balls around direct-hyperbolic unstable periodic orbits (UPOs), then
% chooses switching iterates by minimizing the remaining resolved transfer
% time to the next proxy or final target.
%
% Stable and unstable manifolds are diagnostic overlays.  They are not used
% as a complete explicit polygonal pass network.

clear; close all; clc;
projectFolder = fileparts(mfilename('fullpath'));
restoredefaultpath;
rehash toolboxcache;
addpath(projectFolder);
fprintf('Project folder: %s\n', projectFolder);
fprintf('MATLAB version: %s\n', version);

%% Authoritative implementation and file inventory
% Only this folder is placed on the MATLAB path.  The nested backup/copy
% folder is deliberately not added to the path.  The inventory table below
% is generated from the files that MATLAB can see in this authoritative
% folder.

authoritativeFolder = projectFolder;
authoritativeFunctions = {
    'schroer_ott_default_config'
    'so_standard_map_lifted'
    'so_standard_map_inverse_lifted'
    'so_jacobian'
    'so_jacobian_product'
    'so_enumerate_periodic_orbits'
    'so_construct_first_light_route'
    'so_resolve_connection'
    'so_multistage_targeting'
    'so_replay_stage'
    'so_plot_results'
    };
calledFiles = strings(numel(authoritativeFunctions), 1);
for k = 1:numel(authoritativeFunctions)
    calledFiles(k) = string(which(authoritativeFunctions{k}));
end
calledFileTable = table(string(authoritativeFunctions), calledFiles, ...
    'VariableNames', {'functionName','resolvedPath'})

inventory = local_file_inventory(projectFolder)

%% Mathematical and numerical configuration
% The corrected deterministic first-light case is not the paper Figure-3
% geometry.  It is a nondegenerate development case at k=1.25.  Its source
% centre is outside every retained proxy ball and remains outside those
% balls under every admissible zero-time vertical control.
%
% The route bracket has a declared margin and open boundaries.  This avoids
% accepting or rejecting a resonance merely because a rectangle centre is
% equal to a nominal rotation number to round-off precision.

cfg = schroer_ott_default_config();
cfg.saveFigures = true;
if ~exist(cfg.outputDirectory, 'dir'); mkdir(cfg.outputDirectory); end
if ~exist(cfg.figureDirectory, 'dir'); mkdir(cfg.figureDirectory); end
configurationTable = local_configuration_table(cfg)

%% Standard map, inverse map, Jacobian, and control segment
% The lifted standard map is
%
%   y_{n+1} = y_n - k/(2*pi) sin(2*pi*x_n),
%   x_{n+1} = x_n + y_{n+1}.
%
% The cylinder version wraps only x.  The inverse map is checked by applying
% the forward map followed by the inverse map.  The Jacobian determinant is
% one, so the map is area-preserving in the lifted variables.

zCheck = [0.31; 0.44];
zForward = so_standard_map_lifted(zCheck, cfg);
zBack = so_standard_map_inverse_lifted(zForward, cfg);
J = so_jacobian(zCheck, cfg);
h = 1e-7;
Jfd = [(so_standard_map_lifted(zCheck + [h; 0], cfg) - so_standard_map_lifted(zCheck - [h; 0], cfg)) / (2*h), ...
       (so_standard_map_lifted(zCheck + [0; h], cfg) - so_standard_map_lifted(zCheck - [0; h], cfg)) / (2*h)];
mapValidationTable = table( ...
    norm(zBack - zCheck), ...
    det(J), ...
    norm(J - Jfd), ...
    cfg.controlAmplitude, ...
    2*cfg.controlAmplitude, ...
    'VariableNames', {'forwardInverseError','jacobianDeterminant','jacobianFiniteDifferenceError','uMax','delta'})

%% UPO enumeration
% Periodic chains are solved using fixed-(p,m) Newton equations in the
% lifted map,
%
%   F_lift^p(z) - z - [m; 0] = 0.
%
% Multiple distinct chains with the same (p,m) are retained.  Lower-period
% duplicates are rejected, and cyclic/cylindrical deduplication prevents
% repeated copies of the same chain.  Stability is classified from the
% trace of the monodromy matrix.

catalogue = so_enumerate_periodic_orbits(cfg);
upoCatalogue = catalogue.byOmega
retainedDirectHyperbolic = catalogue.byOmega(catalogue.byOmega.classification == "direct-hyperbolic", :)
diagnosticFractions = catalogue.diagnosticFractions
maximumAverageYError = max(catalogue.byOmega.averageYError)

%% Route construction and source-to-proxy audit
% The route is chosen algorithmically from direct-hyperbolic chains whose
% rotation numbers fall inside the declared source-target y bracket.  No
% chain identifier is hard-coded into the route.
%
% The audit table explicitly checks the corrected nondegeneracy condition:
% the source starts outside all initial proxy balls, and the entire
% admissible zero-time control segment also stays outside them.

route = so_construct_first_light_route(catalogue, cfg);
routeTable = local_route_table(route)
sourceProxyAudit = so_source_proxy_audit(cfg, route);
sourceProxyAudit.byProxy
auditSummary = struct2table(rmfield(sourceProxyAudit, 'byProxy'), 'AsArray', true)

%% Forward-backward connection solver
% A connection from a source point to a target component is found by a
% diagonal search in increasing total time.  For each split (n_forward,
% n_backward), the forward image of the vertical control segment and the
% backward image of the target boundary are adaptively refined and
% intersected.  Intersections are refined inside their original segment
% brackets.  Direct zero-time containment is treated separately and chooses
% the admissible control with minimum absolute value.
%
% The first provisional connection below is computed only to expose the
% connection data before the full multistage solve.

sourceState = so_rectangle_center(cfg.sourceRectangle);
sourceState(1) = so_wrap_x(sourceState(1));
[firstConnection, firstProfile] = so_resolve_connection(sourceState, route.targetComponents{1}, cfg, so_new_performance_profile());
firstConnectionSummary = local_connection_table(firstConnection)
firstProfile

%% Multistage targeting and adaptive switch rule
% The full targeting routine first resolves a provisional connection to the
% current proxy.  It then tests each possible switching index j along that
% planned path.  The selected j minimizes
%
%   J(j) = j + tau_resolved(w_j -> next target).
%
% If j=0 is selected, the current proxy is skipped because the calculated
% route to the next target is no worse than entering the proxy first.  The
% skip is therefore an objective-based decision, not evidence that the
% source began inside the proxy.

targetingTimer = tic;
result = so_multistage_targeting(cfg, catalogue, route);
result.diagnosticManifolds = so_compute_diagnostic_manifolds(route, cfg);
targetingRuntime = toc(targetingTimer);
fprintf('Corrected deterministic targeting runtime: %.3f seconds\n', targetingRuntime);
stageSummary = so_stage_summary_table(result)
switchTable = local_switch_table(result.switchEvaluations)
executedControlTable = result.executedControls

%% Independent replay of each accepted stage
% Each accepted executed segment stores the pre-control state, the
% instantaneous control, the post-control state, the planned path, the
% switch/final state, and pathwise replay errors.  Replay starts from the
% stored pre-control state, applies the stored y-kick once, and advances the
% lifted standard map one step at a time.  It does not call the connection
% solver and does not reconstruct the control.

replayTable = local_replay_table(result)

%% Controlled trajectory figures
% The figures are generated inline.  PNG copies are optional, but the
% important thesis figures are also exported as SVG and EPS in the figure
% directory.  The trajectory figure preserves instantaneous control jumps
% and inserts NaN breaks across periodic x-boundary jumps.

figureHandles = so_plot_results(result);
exportedFigureInventory = local_export_inventory(cfg.figureDirectory)

%% Uncontrolled propagation from the same source point
% The uncontrolled baseline is propagated from the same source centre with
% no controls.  Hitting is checked against the final target rectangle over a
% declared finite observation horizon.

uncontrolledHorizon = 5000;
uncontrolledComparison = local_uncontrolled_comparison(sourceState, cfg, uncontrolledHorizon)

%% Equation (1) diagnostic
% The paper's notation for target size is ambiguous because it refers to a
% target radius and later to a ball diameter.  The code declares the
% convention in cfg.equationOne.targetSizeConvention.  The table also shows
% the predicted value obtained under the alternative convention.  L=1 is
% labelled as a declared order-one phase-space scale, not a fitted length.

equationOneDiagnostics = result.equationOneDiagnostics

%% Corrected validation suite
% The tests include independent negative replay checks: corrupting a stored
% planned path or corrupting a stored control must fail replay.  The test
% names are distinct and no longer call other tests as aliases.

testResults = test_schroer_ott_first_light()
testSummary = groupsummary(testResults, "passed")

%% Three-level numerical-resolution study
% This compact resolution study is project-specific.  It varies only curve
% and intersection-resolution settings for the corrected deterministic
% first-light case.  It does not change k, source, target, control bound,
% route policy, or iteration horizons.

resolutionStudy = so_run_resolution_study(true);
resolutionStudy.comparison

%% Estimated paper Figure-3 comparison case
% The source and target rectangles for the paper comparison are estimated
% from the rendered Figure 3 axes because the paper does not tabulate their
% exact coordinates.  The seven route rotations are selected from the
% direct-hyperbolic chains matching the full-width resonances visible in
% the figure: 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, and 3/4.
%
% The published 125-132 iteration range is reported descriptively.  It is
% not used as an acceptance criterion.

[paperCfg, paperMetadata] = so_paper_comparison_config();
paperMetadata
paperCase = so_run_paper_comparison_case(true);
paperComparisonSummary = local_paper_case_table(paperCase)

%% Fifty-source paper-comparison ensemble
% The following command runs exactly 50 uniformly sampled source points in
% the adopted paper-comparison source rectangle.  It stores all successes,
% failures, unresolved cases, exact sampled states, replay errors, and
% comparison counts relative to the published 125-132 range.
%
% This can be computationally expensive because every trial runs the full
% forward-backward route independently.  It is intentionally not replaced by
% three hand-picked points.

paperEnsembleFile = fullfile(paperCfg.outputDirectory, 'paper_50_source_ensemble.mat');
runPaperEnsembleNow = exist('RUN_FULL_PAPER_ENSEMBLE', 'var') && RUN_FULL_PAPER_ENSEMBLE;
if runPaperEnsembleNow || ~isfile(paperEnsembleFile)
    paperEnsemble = so_run_paper_ensemble(true);
else
    loadedEnsemble = load(paperEnsembleFile, 'ensemble');
    paperEnsemble = loadedEnsemble.ensemble;
end
paperEnsemble.summary
paperEnsemble.trials

%% One mild-noise demonstration
% The deterministic map functions are not changed for noise.  Instead, a
% separate replay wrapper applies the already-computed deterministic
% control schedule and adds one stored Gaussian y-kick sequence with
% sigma = 0.01*(2*u_max) = 6e-5.  No retargeting is performed.

noiseDemo = so_run_mild_noise_demo(result, true);
noiseSummary = local_noise_table(noiseDemo)

%% Final validation summary
% The table below collects the main deterministic checks.  Paper-comparison
% and 50-source ensemble outcomes must be interpreted separately because
% the paper geometry is estimated from a figure, not tabulated.

finalValidationSummary = local_final_validation_table(result, sourceProxyAudit, testResults, resolutionStudy, paperCase, paperEnsemble, noiseDemo)

%% Local helper functions
function inventory = local_file_inventory(projectFolder)
files = dir(fullfile(projectFolder, '*'));
isRegular = ~[files.isdir].';
files = files(isRegular);
inventory = table(string({files.name}).', string({files.folder}).', [files.bytes].', ...
    datetime([files.datenum].', 'ConvertFrom', 'datenum'), ...
    'VariableNames', {'name','folder','bytes','lastModified'});
inventory = sortrows(inventory, 'name');
end

function tbl = local_configuration_table(cfg)
names = ["k";"controlAmplitude";"sourceX";"sourceY";"targetX";"targetY"; ...
    "maxForwardIterations";"maxBackwardIterations";"maxTotalTransferTime"; ...
    "proxyTargetRadius";"routeBracketMargin";"routeLowerBoundary";"routeUpperBoundary"; ...
    "intersectionBackend";"equationOneTargetSizeConvention";"equationOneLConvention"];
values = [
    string(cfg.k)
    string(cfg.controlAmplitude)
    sprintf('[%.6f, %.6f]', cfg.sourceRectangle.xMin, cfg.sourceRectangle.xMax)
    sprintf('[%.6f, %.6f]', cfg.sourceRectangle.yMin, cfg.sourceRectangle.yMax)
    sprintf('[%.6f, %.6f]', cfg.targetRectangle.xMin, cfg.targetRectangle.xMax)
    sprintf('[%.6f, %.6f]', cfg.targetRectangle.yMin, cfg.targetRectangle.yMax)
    string(cfg.maxForwardIterations)
    string(cfg.maxBackwardIterations)
    string(cfg.maxTotalTransferTime)
    string(cfg.proxyTargetRadius)
    string(cfg.routeBracket.margin)
    string(cfg.routeBracket.lowerBoundary)
    string(cfg.routeBracket.upperBoundary)
    string(cfg.intersectionBackend)
    string(cfg.equationOne.targetSizeConvention)
    string(cfg.equationOne.LConvention)
    ];
tbl = table(names, values, 'VariableNames', {'quantity','value'});
end

function tbl = local_route_table(route)
chainID = string({route.chains.id}).';
omega = [route.chains.omega].';
period = [route.chains.period].';
winding = [route.chains.winding].';
classification = string({route.chains.classification}).';
componentCount = cellfun(@numel, route.targetComponents).';
tbl = table(chainID, omega, period, winding, classification, componentCount);
end

function tbl = local_connection_table(conn)
tbl = table(conn.success, conn.nForward, conn.nBackward, conn.totalIterations, ...
    conn.control, conn.intersectionResidual, conn.directContainment, ...
    conn.zeroTimePretest, conn.targetContained, conn.timeMinimumCertified, ...
    conn.selectionCertified, string(conn.targetComponent), ...
    'VariableNames', {'success','nForward','nBackward','totalIterations','controlY', ...
    'intersectionResidual','directContainment','zeroTimePretest','targetContained', ...
    'timeMinimumCertified','selectionCertified','targetComponent'});
end

function tbl = local_switch_table(switchEvaluations)
rows = {};
for i = 1:numel(switchEvaluations)
    selected = switchEvaluations(i).selected;
    rows(end + 1, :) = {switchEvaluations(i).stage, switchEvaluations(i).currentTargetID, ...
        switchEvaluations(i).nextTargetID, selected.j, selected.objectiveJ, ...
        selected.finite, switchEvaluations(i).pruned, numel(switchEvaluations(i).probes)}; %#ok<AGROW>
end
tbl = cell2table(rows, 'VariableNames', {'stage','currentTargetID','nextTargetID', ...
    'selectedJ','selectedObjective','selectedFinite','prunedProbeCount','probeCount'});
end

function tbl = local_replay_table(result)
rows = {};
for i = 1:numel(result.executionSegments)
    seg = result.executionSegments(i);
    rows(end + 1, :) = {i, seg.preControlState(1), seg.preControlState(2), ...
        seg.control, seg.postControlState(1), seg.postControlState(2), ...
        seg.iterations, seg.maxPathwiseReplayError, seg.switchStateReplayError, ...
        seg.finalStateReplayError, seg.targetContained, seg.replayPassed}; %#ok<AGROW>
end
tbl = cell2table(rows, 'VariableNames', {'segment','preControlX','preControlY', ...
    'controlY','postControlX','postControlY','iterations','maxPathwiseReplayError', ...
    'switchStateReplayError','finalStateReplayError','targetContained','replayPassed'});
end

function exported = local_export_inventory(figureDirectory)
files = [dir(fullfile(figureDirectory, '*.svg')); dir(fullfile(figureDirectory, '*.eps')); ...
    dir(fullfile(figureDirectory, '*.png'))];
if isempty(files)
    exported = table();
else
    exported = table(string({files.name}).', [files.bytes].', ...
        'VariableNames', {'name','bytes'});
end
end

function tbl = local_uncontrolled_comparison(z0, cfg, horizon)
target = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
z = z0;
hitIteration = NaN;
for n = 1:horizon
    z = so_standard_map_lifted(z, cfg);
    if target.contains(z, cfg.containmentTolerance)
        hitIteration = n;
        break;
    end
end
tbl = table(horizon, ~isnan(hitIteration), hitIteration, z(1), z(2), ...
    'VariableNames', {'observationHorizon','reached','hitIteration','lastX','lastY'});
end

function tbl = local_paper_case_table(paperCase)
if isempty(paperCase.result)
    tbl = table(string(paperCase.classification), false, Inf, 0, false, ...
        'VariableNames', {'classification','targetReached','totalIterations','numberOfControls','targetContained'});
else
    tbl = table(string(paperCase.classification), paperCase.result.targetReached, ...
        paperCase.result.totalExecutedIterations, paperCase.result.numberOfControls, ...
        paperCase.result.targetContained, ...
        'VariableNames', {'classification','targetReached','totalIterations','numberOfControls','targetContained'});
end
end

function tbl = local_noise_table(noiseDemo)
tbl = table(noiseDemo.sigma, noiseDemo.delta, noiseDemo.seed, noiseDemo.noisySuccess, ...
    noiseDemo.totalIterations, noiseDemo.maximumControl, noiseDemo.storedNoiseReplayError, ...
    noiseDemo.finalState(1), noiseDemo.finalState(2), noiseDemo.finalContainment, ...
    noiseDemo.endpointDisplacement, ...
    'VariableNames', {'sigma','delta','seed','noisySuccess','totalIterations', ...
    'maximumControl','storedNoiseReplayError','finalX','finalY','finalContainment', ...
    'endpointDisplacement'});
end

function tbl = local_final_validation_table(result, audit, testResults, resolutionStudy, paperCase, paperEnsemble, noiseDemo)
quantity = [
    "deterministicTargetReached"
    "sourceOutsideInitialProxies"
    "controlSegmentOutsideInitialProxies"
    "anyPositiveSwitch"
    "allControlsWithinBound"
    "maxReplayError"
    "finalContainment"
    "testSuitePassed"
    "resolutionLevelsSucceeded"
    "paperCaseClassification"
    "paperEnsembleAttempts"
    "paperEnsembleSuccesses"
    "noiseStoredReplayError"
    ];
value = [
    string(result.targetReached)
    string(audit.sourceOutsideAllInitialProxies)
    string(audit.controlSegmentOutsideAllInitialProxies)
    string(any(result.executedControls.stageIterations > 0))
    string(all(result.executedControls.controlWithinBound))
    string(max(result.executedControls.maxPathwiseReplayError))
    string(result.targetContained)
    string(all(testResults.passed))
    string(all(resolutionStudy.comparison.success))
    string(paperCase.classification)
    string(paperEnsemble.summary.attempts)
    string(paperEnsemble.summary.successes)
    string(noiseDemo.storedNoiseReplayError)
    ];
tbl = table(quantity, value);
end
