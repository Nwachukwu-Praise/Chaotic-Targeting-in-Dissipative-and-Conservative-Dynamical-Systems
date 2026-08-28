function report = compare_shinbrot_search_versions(mapData)
%COMPARE_SHINBROT_SEARCH_VERSIONS Three searches, one source, side by side.
%
%   Runs all three implementations from accepted crossing 150 and reports
%   what each returns, what it cost, and how it relates to the published
%   Table 2 values.
%
%     pre-2026-08-09   reconstruction.  Discards a same-signature interval
%                      whose endpoints do not bracket X_t, so it can return
%                      a larger perturbation than necessary, or miss a
%                      horizon entirely.
%     post-2026-08-09  current file.  Fixes that, but adds an unconditional
%                      257-point sweep at every horizon and estimated
%                      Lipschitz pruning, neither of which is in the paper.
%     coverage         keeps the fix, removes the additions, and gates each
%                      horizon on the paper's own condition X_t in DeltaX_n.
%
%   The comparison to look at is not only which p is returned, but the
%   parameter-evaluation count, which is hardware independent and is the
%   quantity Shinbrot et al.'s description of the method as "quick" can
%   reasonably be read against.

if nargin < 1 || isempty(mapData)
    mapData = local_nominal_map();
end

params = local_params();
target = local_target(params);
control = local_control();

sourceIndex = 150;
sourceState = mapData.states(:, sourceIndex);
while abs(sourceState(1) - target.x) <= target.tolerance
    sourceIndex = sourceIndex + 1;
    sourceState = mapData.states(:, sourceIndex);
end

published = struct('p', 0.086719, 'error', 0.0069203, 'horizon', 3, ...
    'searchTime', 1.3707);

fprintf('\n=== Shinbrot search versions, source %d (X_s = %.8g) ===\n', ...
    sourceIndex, sourceState(1));
fprintf('published Table 2: p = %.6g, error = %.6g, n = %d, %.4f s\n\n', ...
    published.p, published.error, published.horizon, published.searchTime);

versions = { ...
    'pre-2026-08-09',  @run_shinbrot_pre20260809_bisection; ...
    'post-2026-08-09', @run_shinbrot_discontinuity_aware_bisection; ...
    'coverage',        @run_shinbrot_coverage_bisection};

rows = {};
searches = cell(size(versions, 1), 1);
for v = 1:size(versions, 1)
    [row, searchOut] = local_run(versions{v, 1}, versions{v, 2}, ...
        sourceState, target, params, control, published);
    rows(end + 1, :) = row; %#ok<AGROW>
    searches{v} = searchOut;
end

report.table = cell2table(rows, 'VariableNames', {'version','found', ...
    'selectedP','absP','finalTargetError','horizon','searchSeconds', ...
    'parameterEvaluations','crossingPropagations','bracketBisections', ...
    'exploratorySubdivisions','replayHit','replayTargetError', ...
    'matchesPublished'});
report.searches = searches;
report.published = published;
report.sourceIndex = sourceIndex;
report.sourceState = sourceState;

disp(report.table);

fprintf('every returned perturbation is replayed independently from the\n');
fprintf('three-dimensional source state before it is accepted.\n\n');

matched = report.table.matchesPublished;
if any(matched)
    fprintf('reproduces Table 2 : %s\n', ...
        strjoin(report.table.version(matched), ', '));
else
    fprintf(['reproduces Table 2 : none.  The reconstruction is then not\n', ...
        'the original; traversal order is the likeliest cause.\n']);
end

evaluations = report.table.parameterEvaluations;
if all(isfinite(evaluations)) && evaluations(1) > 0
    fprintf('\nparameter evaluations relative to pre-2026-08-09:\n');
    for v = 1:height(report.table)
        fprintf('  %-17s %6d  (%.1fx)\n', report.table.version{v}, ...
            evaluations(v), evaluations(v) / evaluations(1));
    end
end

coverageRow = find(strcmp(report.table.version, 'coverage'), 1);
if ~isempty(coverageRow) && ~isempty(searches{coverageRow})
    s = searches{coverageRow};
    if isfield(s, 'horizonRejectedByImageSpan')
        rejected = find(s.horizonRejectedByImageSpan(:).');
        searched = find(s.horizonSearched(:).');
        fprintf('\ncoverage gate, Equation (5):\n');
        fprintf('  probe samples per horizon        : %d\n', s.imageProbeSamples);
        fprintf('  horizons rejected on image span  : %s\n', mat2str(rejected));
        fprintf('  horizons actually searched       : %s\n', mat2str(searched));
        fprintf(['  A rejected horizon is one where the probed image set\n', ...
            '  does not reach X_t, which is the paper''s own condition for\n', ...
            '  advancing.  It is not an exhausted search.\n']);
    end
end
fprintf('=========================================================\n');
end

% -------------------------------------------------------------------------

function [row, search] = local_run(label, searchFunction, sourceState, ...
    target, params, control, published)
search = local_blank_search();
bestApprox = struct('p', NaN, 'horizon', NaN, 'x', NaN, 'error', Inf);
stats = struct('parameterEvaluations', 0, 'crossingPropagations', 0, ...
    'bisectionIterations', 0, 'crossingSignatureEvaluations', 0);

timer = tic;
try
    [search, ~, stats] = searchFunction(search, bestApprox, stats, ...
        sourceState, target, params, control, -control.deltaP, ...
        control.deltaP, control.maxSearchCrossings);
catch err
    fprintf('%s FAILED: %s\n', label, err.message);
    row = {label, false, NaN, NaN, NaN, NaN, toc(timer), NaN, NaN, NaN, ...
        NaN, false, NaN, false};
    return;
end
elapsed = toc(timer);

replay = local_replay(sourceState, target, params, ...
    search.selectedP, search.horizon);

bracketBisections = local_field(search, 'bracketBisectionIterations', ...
    stats.bisectionIterations);
exploratory = local_field(search, 'exploratorySubdivisions', NaN);

matches = isfinite(search.selectedP) && ...
    abs(search.selectedP - published.p) <= 5e-6 && ...
    abs(search.finalTargetError - published.error) <= 5e-7;

row = {label, search.found, search.selectedP, abs(search.selectedP), ...
    search.finalTargetError, search.horizon, elapsed, ...
    stats.parameterEvaluations, stats.crossingPropagations, ...
    bracketBisections, exploratory, replay.hit, replay.targetError, matches};
end

function replay = local_replay(sourceState, target, params, pControl, horizon)
replay.hit = false;
replay.targetError = NaN;
if ~isfinite(pControl) || ~isfinite(horizon) || horizon < 1
    return;
end
currentState = sourceState(:);
for crossing = 1:horizon
    segment = next_valid_section_crossing(currentState, params, pControl);
    if ~segment.success
        return;
    end
    currentState = segment.eventState;
end
replay.targetError = abs(currentState(1) - target.x);
replay.hit = replay.targetError <= target.tolerance;
end

function value = local_field(s, name, defaultValue)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = defaultValue;
end
end

function search = local_blank_search()
search.found = false;
search.selectedP = NaN;
search.horizon = NaN;
search.targetX = NaN;
search.bestDistance = Inf;
search.finalTargetError = Inf;
search.selectionMethod = '';
search.bracketLow = NaN;
search.bracketHigh = NaN;
search.finalBracketBisectionIterations = 0;
search.bisectionIterationsUsed = 0;
search.plotHorizon = NaN;
search.gridEvaluations = 0;
search.failureCode = 'none';
search.failureReason = '';
end

function params = local_params()
params.sigma = 10;
params.rho = 28;
params.beta = 8/3;
params.zSection = 26.921;
params.xSectionMin = 8.0;
params.crossingDirection = +1;
params.odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
params.eventDisableTime = 1e-6;
params.maxStepOffSection = 1e-5;
params.sectionTolerance = 1e-9;
params.maxTimeToNextValidSection = 100;
params.maxRejectedCrossings = 100;
end

function target = local_target(params)
target.x = 13.729;
target.state = [13.729; 19.585; params.zSection];
target.tolerance = 0.008;
end

function control = local_control()
control.deltaP = 0.1;
control.maxSearchCrossings = 12;
control.numParameterSamples = 81;
control.bisectionIterations = 24;
control.maxDiscontinuityIsolationDepth = 24;
control.maxDiscontinuityIntervals = 256;
end

function mapData = local_nominal_map()
cacheFile = fullfile(fileparts(mfilename('fullpath')), 'shinbrot_nominal_map.mat');
if isfile(cacheFile)
    cached = load(cacheFile, 'mapData');
    mapData = cached.mapData;
    return;
end
fprintf('Generating the nominal return map (1200 crossings)...\n');
mapData = generate_return_map([1; 1; 1], 1200, local_params(), 0);
save(cacheFile, 'mapData', '-v7.3');
end
