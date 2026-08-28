function report = verify_pre20260809_reconstruction(mapData)
%VERIFY_PRE20260809_RECONSTRUCTION Does the reconstruction reproduce Table 2?
%
%   The reconstruction in run_shinbrot_pre20260809_bisection is rebuilt from
%   a prose description, not recovered from a backup.  It is only the
%   original if it returns the values the report already publishes:
%
%       p = 0.086719,  final target error = 0.0069203,  horizon n = 3
%
%   from the source at accepted crossing 150.  This script runs both the
%   reconstruction and the current search from that source and prints the
%   comparison.  It does not decide anything; it establishes which of the
%   two searches produced the published table.

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

published.p = 0.086719;
published.error = 0.0069203;
published.horizon = 3;
published.searchTime = 1.3707;

fprintf('\n=== Which search produced Table 2? ===\n');
fprintf('source index %d, X_s = %.8g\n', sourceIndex, sourceState(1));
fprintf('published: p = %.6g, error = %.6g, n = %d, search %.4f s\n\n', ...
    published.p, published.error, published.horizon, published.searchTime);

rows = {};
rows(end + 1, :) = local_run('reconstruction (pre 2026-08-09)', ...
    @run_shinbrot_pre20260809_bisection, sourceState, target, params, control);
rows(end + 1, :) = local_run('current (post 2026-08-09)', ...
    @run_shinbrot_discontinuity_aware_bisection, sourceState, target, params, control);

report.table = cell2table(rows, 'VariableNames', {'implementation','found', ...
    'selectedP','finalTargetError','horizon','searchSeconds', ...
    'parameterEvaluations','pMatchesPublished','errorMatchesPublished'});
report.published = published;
report.sourceState = sourceState;
report.sourceIndex = sourceIndex;

disp(report.table);

matched = report.table.pMatchesPublished & report.table.errorMatchesPublished;
if any(matched)
    fprintf('Table 2 was produced by: %s\n', ...
        report.table.implementation{find(matched, 1)});
else
    fprintf(['Neither implementation reproduces Table 2.\n', ...
        'The reconstruction is then not the original.  The traversal order\n', ...
        'is the likeliest discrepancy: it is not recorded in the surviving\n', ...
        'description, and a successful horizon returns at the first hit, so\n', ...
        'order selects among several admissible perturbations.  Try popping\n', ...
        'the right half first in explore_horizon before concluding further.\n']);
end
fprintf('======================================\n');

    function row = local_run(label, searchFunction, zSource, tgt, prm, ctl)
        search = local_blank_search();
        bestApprox = struct('p', NaN, 'horizon', NaN, 'x', NaN, 'error', Inf);
        stats = struct('parameterEvaluations', 0, 'crossingPropagations', 0, ...
            'bisectionIterations', 0, 'crossingSignatureEvaluations', 0);
        timer = tic;
        [search, ~, stats] = searchFunction(search, bestApprox, stats, ...
            zSource, tgt, prm, ctl, -ctl.deltaP, ctl.deltaP, ...
            ctl.maxSearchCrossings);
        elapsed = toc(timer);
        pMatch = isfinite(search.selectedP) && ...
            abs(search.selectedP - published.p) <= 5e-6;
        eMatch = isfinite(search.finalTargetError) && ...
            abs(search.finalTargetError - published.error) <= 5e-7;
        row = {label, search.found, search.selectedP, search.finalTargetError, ...
            search.horizon, elapsed, stats.parameterEvaluations, pMatch, eMatch};
    end
end

% -------------------------------------------------------------------------

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
