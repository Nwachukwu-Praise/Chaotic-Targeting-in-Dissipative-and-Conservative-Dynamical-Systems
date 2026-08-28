function diagnostics = shinbrot_chapter4_diagnostics(parts, settings)
%SHINBROT_CHAPTER4_DIAGNOSTICS Generate the numbers Section 4 is missing.
%
%   Run from the Shinbrot targeting folder:
%
%       shinbrot_chapter4_diagnostics                      % all parts
%       shinbrot_chapter4_diagnostics(["replay" "noise"])  % a subset
%
%   Everything printed also goes to shinbrot_chapter4_diagnostics.log.
%
%   Parts
%     "map"        nominal return map: mean return time, mean log slope of
%                  equation (27), and the one-dimensionality thickness.
%                  Regenerating the 1200-crossing map is the slow step; it is
%                  cached in shinbrot_nominal_map.mat and reused.
%     "resolution" the declared search resolution of equation (12), the
%                  horizon beyond which it no longer covers the expected
%                  feature width, and the discontinuity counts for the
%                  principal demonstration.
%     "replay"     replay independence and the two negative controls.
%     "noise"      Tbar, N, sigma_*, the delta conversion of equation (23),
%                  the adaptive versus fixed-step discrepancy, and Table 4
%                  with Wilson intervals and the retargeting-failure split.
%
%   This produces no new physics.  It reports quantities Subsections 3.1.5
%   and 3.1.9 already promise, from results that already exist.

if nargin < 1 || isempty(parts)
    parts = ["map", "resolution", "replay", "noise"];
end
parts = string(parts);
if nargin < 2
    settings = struct();
end

here = fileparts(mfilename('fullpath'));
addpath(here);
cd(here);

logFile = fullfile(here, 'shinbrot_chapter4_diagnostics.log');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('===============================================================\n');
fprintf('Shinbrot chapter 4 diagnostics\n');
fprintf('date           : %s\n', string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf('MATLAB version : %s\n', version);
fprintf('folder         : %s\n', here);
fprintf('===============================================================\n');

% ------------------------------------------------------------ configuration
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

target.x = 13.729;
target.state = [13.729; 19.585; params.zSection];
target.tolerance = 0.008;

control.deltaP = 0.1;
control.maxSearchCrossings = 12;
control.numParameterSamples = 81;
control.bisectionIterations = 24;
control.maxDiscontinuityIsolationDepth = 24;
control.maxDiscontinuityIntervals = 256;
control.minCertifiedDepth = local_option(settings, 'minCertifiedDepth', 8);
control.searchMethod = verified_bisection_identifier();

noiseDt = local_option(settings, 'noiseDt', 1e-3);
sigmaNoiseValues = local_option(settings, 'sigmaNoiseValues', ...
    [0, 0.01, 0.05, 0.1, 0.5, 1.0]);
numMapCrossings = local_option(settings, 'numMapCrossings', 1200);
sourceIndex = local_option(settings, 'sourceIndex', 150);
noiseTrialsCsv = local_option(settings, 'noiseTrialsCsv', ...
    'noise_trials_bisection_section_retarget.csv');
mapCacheFile = fullfile(here, 'shinbrot_nominal_map.mat');

diagnostics = struct();
diagnostics.params = params;
diagnostics.target = target;
diagnostics.control = control;

% ---------------------------------------------------------- nominal map
needMap = any(ismember(parts, ["map", "resolution", "replay", "noise"]));
mapData = [];
if needMap
    if isfile(mapCacheFile)
        fprintf('\nLoading cached nominal map: %s\n', mapCacheFile);
        cached = load(mapCacheFile, 'mapData', 'numMapCrossings');
        if isfield(cached, 'mapData') && ...
                size(cached.mapData.states, 2) >= numMapCrossings
            mapData = cached.mapData;
        end
    end
    if isempty(mapData)
        fprintf('\nGenerating the nominal return map (%d crossings)...\n', ...
            numMapCrossings);
        tMap = tic;
        mapData = generate_return_map([1; 1; 1], numMapCrossings, params, 0);
        fprintf('  done in %.1f s, %d crossings\n', toc(tMap), ...
            size(mapData.states, 2));
        save(mapCacheFile, 'mapData', 'numMapCrossings', 'params', '-v7.3');
    end
    diagnostics.mapData = mapData;
end

sourceState = mapData.states(:, sourceIndex);
while abs(sourceState(1) - target.x) <= target.tolerance
    sourceIndex = sourceIndex + 1;
    sourceState = mapData.states(:, sourceIndex);
end
diagnostics.sourceIndex = sourceIndex;
diagnostics.sourceState = sourceState;
fprintf('source index %d, X_s = %.8g\n', sourceIndex, sourceState(1));

% ------------------------------------------------------------- part: map
if any(parts == "map")
    fprintf('\n\n############ PART: nominal return map ############\n');
    try
        diagnostics.returnMap = shinbrot_return_map_diagnostics( ...
            mapData, target, struct('burnIn', 100));
    catch err
        fprintf('FAILED: %s\n%s\n', err.message, ...
            getReport(err, 'extended', 'hyperlinks', 'off'));
    end
end

% ------------------------------------------------------ part: resolution
if any(parts == "resolution")
    fprintf('\n\n############ PART: search resolution ############\n');
    try
        diagnostics.resolution = local_resolution_report( ...
            sourceState, target, params, control, diagnostics);
    catch err
        fprintf('FAILED: %s\n%s\n', err.message, ...
            getReport(err, 'extended', 'hyperlinks', 'off'));
    end
end

% ---------------------------------------------------------- part: replay
if any(parts == "replay")
    fprintf('\n\n############ PART: replay verification ############\n');
    try
        diagnostics.replay = test_shinbrot_replay_negative_control( ...
            sourceState, target, params, control);
    catch err
        fprintf('FAILED: %s\n%s\n', err.message, ...
            getReport(err, 'extended', 'hyperlinks', 'off'));
    end
end

% ----------------------------------------------------------- part: noise
if any(parts == "noise")
    fprintf('\n\n############ PART: noise calibration ############\n');
    try
        adaptiveError = NaN;
        searchNow = search_parameter_to_target(sourceState, target, params, ...
            control, verified_bisection_identifier());
        if searchNow.found
            adaptiveError = searchNow.finalTargetError;
        end

        fixedStepError = NaN;
        if isfile(noiseTrialsCsv)
            trials = readtable(noiseTrialsCsv);
            zeroMask = trials.SigmaNoise == 0;
            if any(zeroMask)
                fixedStepError = mean(trials.FinalTargetError(zeroMask), 'omitnan');
            end
        end

        calibration = shinbrot_noise_calibration(mapData, target, noiseDt, ...
            sigmaNoiseValues, struct('burnIn', 100, ...
            'adaptiveFinalError', adaptiveError, ...
            'fixedStepFinalError', fixedStepError));
        diagnostics.calibration = calibration;

        if isfile(noiseTrialsCsv)
            diagnostics.noiseSummary = shinbrot_noise_summary_extended( ...
                noiseTrialsCsv, calibration);
        else
            fprintf('\n%s not found; extended Table 4 skipped.\n', noiseTrialsCsv);
        end
    catch err
        fprintf('FAILED: %s\n%s\n', err.message, ...
            getReport(err, 'extended', 'hyperlinks', 'off'));
    end
end

outFile = fullfile(here, 'shinbrot_chapter4_diagnostics.mat');
save(outFile, 'diagnostics', '-v7.3');
fprintf('\n===============================================================\n');
fprintf('saved: %s\n', outFile);
fprintf('log  : %s\n', logFile);
fprintf('===============================================================\n');
diary off;
end

% =========================================================================

function report = local_resolution_report(sourceState, target, params, ...
    control, diagnostics)
%LOCAL_RESOLUTION_REPORT Equation (12) resolution and where it stops covering.
%
%   Equation (12) gives the width of the narrowest continuous feature of
%   X_n(p) the search is guaranteed to have sampled.  Subsection 3.1.5 notes
%   that features contract roughly as Delta p * exp(-lambda n), so the
%   guarantee expires at a finite horizon.  Reporting the resolution without
%   that horizon leaves "exhausted at the declared resolution" untestable.

d = control.minCertifiedDepth;
fullWidth = 2 * control.deltaP;
resolution = fullWidth / 2^d;

report.minCertifiedDepth = d;
report.admissibleIntervalWidth = fullWidth;
report.certifiedParameterResolution = resolution;

if isfield(diagnostics, 'returnMap') && ...
        isfield(diagnostics.returnMap, 'meanLogSlope') && ...
        isfinite(diagnostics.returnMap.meanLogSlope)
    lambda = diagnostics.returnMap.meanLogSlope;
    lambdaSource = 'measured <ln|F''|> from the nominal map';
else
    lambda = 0.9;
    lambdaSource = 'assumed 0.9 (run the "map" part for a measured value)';
end
report.lambda = lambda;
report.lambdaSource = lambdaSource;
report.crossoverHorizon = log(fullWidth / resolution) / lambda;

search = search_parameter_to_target(sourceState, target, params, control, ...
    verified_bisection_identifier());
report.search = search;

fprintf('declared resolution, equation (12)\n');
fprintf('  admissible interval width 2*dp   : %.6g\n', fullWidth);
fprintf('  minCertifiedDepth d              : %d\n', d);
fprintf('  dp_res = 2*dp / 2^d              : %.6g\n', resolution);
if isfield(search, 'certifiedParameterResolution')
    fprintf('  reported by the search           : %.6g\n', ...
        search.certifiedParameterResolution);
end
fprintf('  lambda used                      : %.4f  (%s)\n', lambda, lambdaSource);
fprintf('  feature width ~ 2*dp*exp(-lambda n) falls below dp_res at n = %.2f\n', ...
    report.crossoverHorizon);
fprintf(['  Beyond that horizon the search is in its adaptive-descent\n', ...
    '  phase: credible at the declared resolution, not certified.\n']);

fprintf('\nprincipal demonstration, discontinuity awareness\n');
fprintf('  found                            : %d\n', search.found);
fprintf('  selected p                       : %.8g\n', search.selectedP);
fprintf('  horizon                          : %d\n', search.horizon);
fprintf('  discontinuities detected         : %s\n', ...
    local_show(search, 'discontinuitiesDetected'));
fprintf('  discontinuity reductions         : %s\n', ...
    local_show(search, 'discontinuityReductions'));
fprintf('  continuous intervals examined    : %s\n', ...
    local_show(search, 'continuousIntervalsExamined'));
fprintf('  crossing signature evaluations   : %s\n', ...
    local_show(search, 'crossingSignatureEvaluations'));
fprintf('  continuous intervals subdivided  : %s\n', ...
    local_show(search, 'continuousIntervalsSubdivided'));
fprintf('  resolution floor reached         : %s\n', ...
    local_show(search, 'parameterResolutionFloorReached'));
fprintf('  interval budget exhausted        : %s\n', ...
    local_show(search, 'intervalBudgetExhausted'));
fprintf(['\n  A nonzero discontinuity count is the evidence that the\n', ...
    '  discontinuity awareness did work.  A zero count on this source\n', ...
    '  would mean the distinguishing feature of the implementation was\n', ...
    '  never exercised in the principal demonstration, which is worth\n', ...
    '  knowing either way.\n']);
end

function s = local_show(structure, fieldName)
if isfield(structure, fieldName)
    value = structure.(fieldName);
    if islogical(value)
        s = sprintf('%d', value);
    elseif isnumeric(value) && isscalar(value)
        s = sprintf('%g', value);
    else
        s = '(non-scalar)';
    end
else
    s = 'not recorded by this search path';
end
end

function value = local_option(settings, name, defaultValue)
if isstruct(settings) && isfield(settings, name) && ~isempty(settings.(name))
    value = settings.(name);
else
    value = defaultValue;
end
end
