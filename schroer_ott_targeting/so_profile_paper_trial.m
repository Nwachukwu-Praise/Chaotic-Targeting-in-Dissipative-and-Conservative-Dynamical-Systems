function diagnostics = so_profile_paper_trial(parts, trial, runtimeLimitSeconds)
%SO_PROFILE_PAPER_TRIAL Diagnose the paper-geometry ensemble timeout.
%
%   diagnostics = SO_PROFILE_PAPER_TRIAL(parts, trial, runtimeLimitSeconds)
%
%   The initial 50-source paper-comparison ensemble returned 0/50 within
%   its declared budget.  The saved rows classify 45 sources as
%   'runtime_limit_exceeded', 3 as provisional connection failures, and 2 as
%   exhausted switch-probe searches.  This script diagnoses ONE source
%   state before any scientific parameter is changed.
%
%   parts  : any subset of ["census" "hotspots" "stages"], default all.
%              census   - direct cost of building each forward/backward
%                         curve of the paper geometry, per iterate count.
%                         Fast.  Tests the hypothesis that backward
%                         preimages of a proxy circle blow up with n.
%              hotspots - MATLAB code profiler over one capped trial, to
%                         rank functions by self and total time.
%              stages   - one full trial with a generous cap, reporting
%                         per-stage runtime and how far the route gets.
%   trial  : which of the 50 stored ensemble sources to use, default 1.
%            The same seed (1729) and sampling rule as so_run_paper_ensemble
%            are used, so the source state matches the ensemble exactly.
%   runtimeLimitSeconds : cap for the "stages" part, default 1800.
%
%   Nothing here changes the targeting mathematics.  It only measures.

if nargin < 1 || isempty(parts)
    parts = ["census", "hotspots", "stages"];
end
parts = string(parts);
if nargin < 2 || isempty(trial)
    trial = 1;
end
if nargin < 3 || isempty(runtimeLimitSeconds)
    runtimeLimitSeconds = 1800;
end

[cfg, paper] = so_paper_comparison_config();
cfg.verbose = false;
cfg.saveFigures = false;

rng(1729);
sourceStates = [paper.source.xMin + (paper.source.xMax - paper.source.xMin) .* rand(1, 50); ...
    paper.source.yMin + (paper.source.yMax - paper.source.yMin) .* rand(1, 50)];
z0 = sourceStates(:, trial);

fprintf('\n================ paper-trial profile ================\n');
fprintf('trial              : %d\n', trial);
fprintf('source state       : [%.15g; %.15g]\n', z0(1), z0(2));
fprintf('k                  : %g\n', cfg.k);
fprintf('control amplitude  : %g\n', cfg.controlAmplitude);
fprintf('maxForward         : %d\n', cfg.maxForwardIterations);
fprintf('maxBackward        : %d\n', cfg.maxBackwardIterations);
fprintf('maxTotalTransfer   : %d\n', cfg.maxTotalTransferTime);
fprintf('proxy radius       : %g\n', cfg.proxyTargetRadius);
fprintf('curve maxPoints    : %g\n', cfg.curve.maxPoints);
fprintf('curve maxDepth     : %d\n', cfg.curve.maxSubdivisionDepth);

catalogue = so_enumerate_periodic_orbits(cfg);
route = so_construct_omega_route(catalogue, cfg, paper.expectedRotationNumbers, ...
    'estimated Figure-3 seven-resonance route');

fprintf('route rotations    : %s\n', mat2str(route.rotationNumbers, 6));
fprintf('route chain periods: %s\n', mat2str([route.chains.period]));
fprintf('proxy components   : %d (one per phase point, summed over the route)\n', ...
    sum(cellfun(@numel, route.targetComponents)));
fprintf('=====================================================\n');

diagnostics = struct();
diagnostics.trial = trial;
diagnostics.sourceState = z0;
diagnostics.configuration = cfg;
diagnostics.paperMetadata = paper;
diagnostics.route = route;

if any(parts == "census")
    diagnostics.curveCensus = curve_cost_census(z0, route, cfg);
end
if any(parts == "hotspots")
    diagnostics.hotspots = hotspot_profile(z0, catalogue, route, cfg);
end
if any(parts == "stages")
    diagnostics.stageReport = stage_report(z0, catalogue, route, cfg, runtimeLimitSeconds);
end

outFile = fullfile(cfg.outputDirectory, sprintf('paper_trial_profile_%02d.mat', trial));
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
save(outFile, 'diagnostics', '-v7.3');
fprintf('\nsaved: %s\n', outFile);
end

% =========================================================================

function census = curve_cost_census(z0, route, cfg)
%CURVE_COST_CENSUS Cost of each curve the solver must build, per iterate.
fprintf('\n--- part 1: curve cost census ---\n');
fprintf('One row per curve build.  "unresolved" rows are the ones that\n');
fprintf('have hit a budget: they cost the most and resolve nothing.\n\n');

rows = {};

sourceComponent = so_make_control_component(z0, cfg);
for n = 0:cfg.maxForwardIterations
    rows(end + 1, :) = time_one_curve('forward-control-segment', 'control-segment', ...
        sourceComponent, 1, n, cfg); %#ok<AGROW>
end

for r = 1:numel(route.chains)
    comps = route.targetComponents{r};
    comp = comps{1};
    label = sprintf('backward-proxy-omega-%.4g', route.chains(r).omega);
    for n = 0:cfg.maxBackwardIterations
        rows(end + 1, :) = time_one_curve(label, comp.id, comp, -1, n, cfg); %#ok<AGROW>
    end
end

finalComp = route.finalTargetComponents{1};
for n = 0:cfg.maxBackwardIterations
    rows(end + 1, :) = time_one_curve('backward-final-rectangle', finalComp.id, ...
        finalComp, -1, n, cfg); %#ok<AGROW>
end

census = cell2table(rows, 'VariableNames', {'family','componentID','direction', ...
    'iterateCount','pointCount','subdivisionDepth','resolutionStatus', ...
    'maximumGap','maximumMidpointDeviation','curvePointEvaluations','seconds'});

disp(census);

fprintf('total census time            : %.2f s\n', sum(census.seconds));
fprintf('slowest single curve         : %.2f s (%s, n=%d, %d points, %s)\n', ...
    max(census.seconds), ...
    census.family{find(census.seconds == max(census.seconds), 1)}, ...
    census.iterateCount(find(census.seconds == max(census.seconds), 1)), ...
    census.pointCount(find(census.seconds == max(census.seconds), 1)), ...
    census.resolutionStatus{find(census.seconds == max(census.seconds), 1)});
unresolvedMask = ~strcmp(census.resolutionStatus, 'resolved');
fprintf('unresolved curves            : %d of %d\n', sum(unresolvedMask), height(census));
if any(unresolvedMask)
    fprintf('unresolved time              : %.2f s (%.0f%% of census)\n', ...
        sum(census.seconds(unresolvedMask)), ...
        100 * sum(census.seconds(unresolvedMask)) / max(sum(census.seconds), eps));
    fprintf('smallest unresolved iterate  : n=%d\n', min(census.iterateCount(unresolvedMask)));
end

perStage = sum(census.seconds(strcmp(census.family, 'forward-control-segment')));
fprintf('\none forward family build     : %.2f s\n', perStage);
fprintf('switch probes per stage      : up to %d, each rebuilding that family\n', ...
    cfg.maxTotalTransferTime + 1);
fprintf('implied forward cost / stage : %.0f s if every probe is evaluated\n', ...
    perStage * (cfg.maxTotalTransferTime + 1));
fprintf('route stages                 : %d proxies + 1 final\n', numel(route.chains));
end

function row = time_one_curve(family, componentID, component, direction, iterateCount, cfg)
t = tic;
[curve, stats] = so_build_curve(component, direction, iterateCount, cfg);
seconds = toc(t);
row = {family, char(componentID), direction, iterateCount, curve.pointCount, ...
    curve.subdivisionDepth, char(curve.resolutionStatus), curve.maximumGap, ...
    curve.maximumMidpointDeviation, stats.curvePointEvaluations, seconds};
end

% =========================================================================

function hotspots = hotspot_profile(z0, catalogue, route, cfg)
%HOTSPOT_PROFILE MATLAB code profiler over one capped trial.
fprintf('\n--- part 2: function hotspots over one 120 s trial ---\n');

runCfg = cfg;
runCfg.sourceRectangle = struct('xMin', z0(1), 'xMax', z0(1), ...
    'yMin', z0(2), 'yMax', z0(2), 'id', 'paper-source-profile');
timer = tic;
runCfg.runtime.startTime = timer;
runCfg.runtime.limitSeconds = 120;

profile('off');
profile('clear');
profile('on');
result = so_multistage_targeting(runCfg, catalogue, route);
profile('off');
info = profile('info');

entries = info.FunctionTable;
totalTime = arrayfun(@(e) e.TotalTime, entries);
selfTime = arrayfun(@(e) local_self_time(e), entries);
calls = arrayfun(@(e) e.NumCalls, entries);
names = string(arrayfun(@(e) e.FunctionName, entries, 'UniformOutput', false));

[~, order] = sort(totalTime, 'descend');
order = order(1:min(20, numel(order)));

hotspots.table = table(names(order), totalTime(order), selfTime(order), calls(order), ...
    'VariableNames', {'function','totalSeconds','selfSeconds','calls'});
hotspots.result = result;
hotspots.wallClockSeconds = toc(timer);

disp(hotspots.table);
fprintf('wall clock                   : %.2f s\n', hotspots.wallClockSeconds);
fprintf('target reached               : %d\n', local_reached(result));
fprintf('failure category             : %s\n', local_failure(result));
fprintf('stages attempted             : %d\n', numel(result.stagePlans));
end

function t = local_self_time(entry)
t = entry.TotalTime;
if isfield(entry, 'Children') && ~isempty(entry.Children)
    t = t - sum([entry.Children.TotalTime]);
end
t = max(t, 0);
end

% =========================================================================

function report = stage_report(z0, catalogue, route, cfg, runtimeLimitSeconds)
%STAGE_REPORT One trial with a generous cap, per-stage runtime.
fprintf('\n--- part 3: one trial with a %g s cap ---\n', runtimeLimitSeconds);

runCfg = cfg;
runCfg.sourceRectangle = struct('xMin', z0(1), 'xMax', z0(1), ...
    'yMin', z0(2), 'yMax', z0(2), 'id', 'paper-source-stage');
timer = tic;
runCfg.runtime.startTime = timer;
runCfg.runtime.limitSeconds = runtimeLimitSeconds;

result = so_multistage_targeting(runCfg, catalogue, route);
report.wallClockSeconds = toc(timer);
report.result = result;
report.targetReached = local_reached(result);
report.failureCategory = local_failure(result);

rows = {};
for i = 1:numel(result.stagePlans)
    s = result.stagePlans(i);
    conn = s.provisionalConnection;
    rows(end + 1, :) = {i, string(s.targetKind), string(s.targetID), s.rotationNumber, ...
        conn.success, conn.nForward, conn.nBackward, conn.tauResolved, ...
        s.selectedSwitchIndex, s.skipped, s.runtimeSeconds}; %#ok<AGROW>
end
if isempty(rows)
    report.stages = table();
else
    report.stages = cell2table(rows, 'VariableNames', {'stage','targetKind','targetID', ...
        'rotationNumber','connectionSuccess','nForward','nBackward','tauResolved', ...
        'selectedSwitchIndex','skipped','runtimeSeconds'});
    disp(report.stages);
end

p = result.performanceProfile;
fprintf('forward family builds        : %d\n', p.forwardFamilyBuilds);
fprintf('backward family builds       : %d\n', p.backwardFamilyBuilds);
fprintf('backward cache reuses        : %d\n', p.backwardFamilyCacheReuses);
fprintf('refinement calls             : %d\n', p.refinementCalls);
fprintf('unresolved splits            : %d\n', p.unresolvedSplits);
fprintf('switch probes pruned         : %d\n', p.jProbesPruned);
fprintf('peak curve point count       : %d\n', p.peakCurvePointCount);
fprintf('wall clock                   : %.2f s\n', report.wallClockSeconds);
fprintf('target reached               : %d\n', report.targetReached);
fprintf('failure category             : %s\n', report.failureCategory);
fprintf('stages completed             : %d of %d\n', ...
    numel(result.stagePlans), numel(route.chains) + 1);
end

% =========================================================================

function tf = local_reached(result)
tf = isfield(result, 'targetReached') && result.targetReached;
end

function s = local_failure(result)
if isfield(result, 'failureDiagnostics') && ~isempty(result.failureDiagnostics)
    s = string(result.failureDiagnostics(end).failureCategory);
else
    s = "none";
end
end
