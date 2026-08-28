function results = shinbrot_control_bound_sweep(mapData, options)
%SHINBROT_CONTROL_BOUND_SWEEP Subsection 3.1.10, scaling in the control bound.
%
%   results = SHINBROT_CONTROL_BOUND_SWEEP(mapData, options)
%
%   Tests the only quantitative theoretical prediction of the original
%   paper,
%
%       n_t ~ (1/lambda) * ln(1/dp),                              (25)
%
%   by varying the maximum admissible perturbation over the decades
%   dp in [1e-2, 1], directing trajectories to the target from a fixed set
%   of randomly chosen source states at each value, regressing the mean
%   number of accepted section crossings on log10(dp), and converting the
%   fitted slope s to an exponent through
%
%       lambda = ln(10) / |s|.                                    (26)
%
%   Three quantities are compared: the fitted slope s, the published fit
%   s = -2.49 +- 0.09, and the exponent measured directly from the nominal
%   return map as <ln|F'|> of Equation (27), supplied by
%   shinbrot_return_map_diagnostics.
%
%   Subsection 3.1.10 pre-registers three distinguishable outcomes, and the
%   classification below reports which one obtains rather than reconciling
%   the result by adjusting the fit:
%
%     publishedFitReproduced   the fitted slope agrees with the published
%                              value, so the effective exponent in (25) is
%                              not the pointwise mean slope and the
%                              alternative reading is supported;
%     returnMapSlopeReproduced the fitted slope instead agrees with (26)
%                              evaluated at <ln|F'|>, so the present return
%                              map differs from the one underlying the
%                              published fit;
%     neitherReproduced        neither, in which case the discrepancy is
%                              reported and investigated.
%
%   INVESTIGATION OF THE THIRD OUTCOME.  One mechanism is measured here
%   directly.  A map iterate is an ACCEPTED crossing, and the x > 8.0
%   half-plane rejects some plane crossings, so a single map step can span
%   more than one lobe circulation.  That inflates an exponent measured per
%   accepted crossing relative to one measured per plane crossing.  The
%   ratio of plane crossings to accepted crossings is measured on the
%   nominal trajectory and reported, together with lambda rescaled by it.
%
%   HORIZON CAP.  Subsection 3.1.4 fixes H = 12 for the principal case and
%   notes that a cap which binds is indistinguishable in the reported output
%   from an absence of admissible perturbations.  The horizon grows as dp
%   falls, so the cap is raised here and every trial that terminates at the
%   cap is counted, so that a binding cap is visible rather than silent.
%
%   options
%     .deltaPValues    default logspace(-2, 0, 9)
%     .numSources      default 25, following the published protocol
%     .horizonCap      default 25
%     .randomSeed      default 20260814
%     .sourceSeparation default 1.5, matching Equation (9)
%     .candidateRange  default 100:300, matching the ensemble benchmark
%     .returnMapLambda supplied to skip recomputing <ln|F'|>
%     .verbose         default true
%     .saveMatFile     default 'shinbrot_control_bound_sweep_paper_bisection.mat'

if nargin < 1 || isempty(mapData)
    mapData = local_nominal_map();
end
if nargin < 2
    options = struct();
end

deltaPValues   = local_opt(options, 'deltaPValues', logspace(-2, 0, 9));
numSources     = local_opt(options, 'numSources', 25);
horizonCap     = local_opt(options, 'horizonCap', 25);
randomSeed     = local_opt(options, 'randomSeed', 20260814);
separation     = local_opt(options, 'sourceSeparation', 1.5);
candidateRange = local_opt(options, 'candidateRange', 100:300);
verbose        = local_opt(options, 'verbose', true);
saveMatFile    = local_opt(options, 'saveMatFile', ...
    'shinbrot_control_bound_sweep_paper_bisection.mat');

params = local_params();
target = local_target(params);
baseControl = local_control();
baseControl.maxSearchCrossings = horizonCap;

% ---- source states, held fixed across dp -------------------------------
% The same sources are used at every dp so that the regression is paired.
% Drawing fresh sources at each level would add between-source variance to
% the slope, which is the quantity being measured.
eligible = candidateRange(:).';
eligible = eligible(eligible <= size(mapData.states, 2));
eligible = eligible(abs(mapData.x(eligible) - target.x) >= separation);
if numel(eligible) < numSources
    error('shinbrot_control_bound_sweep:TooFewSources', ...
        'Only %d eligible sources; %d requested.', numel(eligible), numSources);
end
rng(randomSeed);
sourceIndices = sort(eligible(randperm(numel(eligible), numSources)));

if verbose
    fprintf('\n=== Scaling in the control bound (Subsection 3.1.10) ===\n');
    fprintf('dp values      : %s\n', mat2str(deltaPValues, 4));
    fprintf('sources        : %d, indices %d..%d, seed %d\n', ...
        numSources, min(sourceIndices), max(sourceIndices), randomSeed);
    fprintf('separation     : |X_s - X_t| >= %g\n', separation);
    fprintf('horizon cap    : %d (raised from the principal case)\n', horizonCap);
    fprintf('total searches : %d\n\n', numel(deltaPValues) * numSources);
end

% ---- sweep --------------------------------------------------------------
trialRows = {};
levelRows = {};
sweepTimer = tic;

for levelIndex = 1:numel(deltaPValues)
    deltaP = deltaPValues(levelIndex);
    control = baseControl;
    control.deltaP = deltaP;

    horizons = NaN(numSources, 1);
    succeeded = false(numSources, 1);
    atCap = false(numSources, 1);
    levelTimer = tic;

    for s = 1:numSources
        sourceIndex = sourceIndices(s);
        sourceState = mapData.states(:, sourceIndex);

        % Follow the dispatcher rather than naming an implementation here.
        % A literal identifier silently pinned this sweep to whichever
        % implementation happened to carry that name, which is how it came
        % to run a search the rest of the project no longer uses.
        search = search_parameter_to_target(sourceState, target, params, ...
            control, verified_bisection_identifier());

        succeeded(s) = search.found;
        if search.found
            horizons(s) = search.horizon;
            atCap(s) = (search.horizon == horizonCap);
        end

        trialRows(end + 1, :) = {deltaP, sourceIndex, mapData.x(sourceIndex), ...
            search.found, search.horizon, search.selectedP, ...
            abs(search.selectedP), search.finalTargetError, ...
            search.parameterEvaluations, search.searchTime, ...
            search.method}; %#ok<AGROW>
    end

    successCount = sum(succeeded);
    if successCount > 0
        successfulHorizons = horizons(succeeded);
        meanHorizon = mean(successfulHorizons);
        medianHorizon = median(successfulHorizons);
        stdHorizon = std(successfulHorizons);
        minHorizon = min(successfulHorizons);
        maxHorizon = max(successfulHorizons);
    else
        meanHorizon = NaN;
        medianHorizon = NaN;
        stdHorizon = NaN;
        minHorizon = NaN;
        maxHorizon = NaN;
    end
    levelRows(end + 1, :) = {deltaP, log10(deltaP), numSources, successCount, ...
        successCount / numSources, meanHorizon, medianHorizon, ...
        stdHorizon, minHorizon, maxHorizon, sum(atCap), ...
        toc(levelTimer)}; %#ok<AGROW>

    if verbose
        fprintf('dp = %8.4g : %2d/%2d succeeded, mean n = %6.3f, at cap %d, %6.1f s\n', ...
            deltaP, successCount, numSources, meanHorizon, sum(atCap), toc(levelTimer));
    end
end

results.trials = cell2table(trialRows, 'VariableNames', {'deltaP','sourceIndex', ...
    'sourceX','found','horizon','selectedP','absP','finalTargetError', ...
    'parameterEvaluations','searchSeconds','searchMethod'});
results.levels = cell2table(levelRows, 'VariableNames', {'deltaP','log10DeltaP', ...
    'sources','successes','successFraction','meanHorizon','medianHorizon', ...
    'stdHorizon','minHorizon','maxHorizon','trialsAtHorizonCap','seconds'});

% ---- regression ---------------------------------------------------------
usable = results.levels.successes > 0 & isfinite(results.levels.meanHorizon);
x = results.levels.log10DeltaP(usable);
y = results.levels.meanHorizon(usable);
fit = local_linear_fit(x, y);

results.fit = fit;
results.fit.lambdaFromSlope = log(10) / abs(fit.slope);
results.fit.lambdaLow  = log(10) / (abs(fit.slope) + 2 * fit.slopeStandardError);
results.fit.lambdaHigh = log(10) / max(abs(fit.slope) - 2 * fit.slopeStandardError, eps);

published.slope = -2.49;
published.slopeUncertainty = 0.09;
published.lambdaRange = log(10) ./ [2.49 + 0.09, 2.49 - 0.09];
results.published = published;

% ---- exponent measured on the return map, Equation (27) -----------------
returnMapLambda = local_opt(options, 'returnMapLambda', []);
if isempty(returnMapLambda)
    rmDiag = shinbrot_return_map_diagnostics(mapData, target, ...
        struct('verbose', false));
    returnMapLambda = rmDiag.meanLogSlope;
    results.returnMapDiagnostics = rmDiag;
end
results.returnMapLambda = returnMapLambda;
results.returnMapImpliedSlope = -log(10) / returnMapLambda;

% ---- plane-to-accepted crossing ratio -----------------------------------
results.crossingRatio = local_crossing_ratio(mapData, params, ...
    sourceIndices(1), 200);
results.lambdaPerPlaneCrossing = returnMapLambda / results.crossingRatio.ratio;

% ---- classification -----------------------------------------------------
slopeLow  = fit.slope - 2 * fit.slopeStandardError;
slopeHigh = fit.slope + 2 * fit.slopeStandardError;
pubLow  = published.slope - published.slopeUncertainty;
pubHigh = published.slope + published.slopeUncertainty;

agreesWithPublished = (slopeLow <= pubHigh) && (pubLow <= slopeHigh);
agreesWithReturnMap = (slopeLow <= results.returnMapImpliedSlope) && ...
    (results.returnMapImpliedSlope <= slopeHigh);

if agreesWithPublished && ~agreesWithReturnMap
    results.classification = "publishedFitReproduced";
elseif agreesWithReturnMap && ~agreesWithPublished
    results.classification = "returnMapSlopeReproduced";
elseif agreesWithPublished && agreesWithReturnMap
    results.classification = "notDiscriminated";
else
    results.classification = "neitherReproduced";
end

results.sourceIndices = sourceIndices;
results.randomSeed = randomSeed;
results.horizonCap = horizonCap;
results.verifiedBisectionIdentifier = verified_bisection_identifier();
results.bisectionImplementation = verified_bisection_implementation();
results.configuration = control_bound_configuration(mapData, params, ...
    target, baseControl, deltaPValues, sourceIndices, randomSeed, ...
    horizonCap, candidateRange, separation, results.returnMapLambda);
results.capBound = any(results.levels.trialsAtHorizonCap > 0);
results.anyLevelIncomplete = any(results.levels.successes < results.levels.sources);
results.runtimeSeconds = toc(sweepTimer);

if ~isempty(saveMatFile)
    save(saveMatFile, 'results', '-v7.3');
end

if ~verbose
    return;
end

fprintf('\n--- levels ---\n');
disp(results.levels);

fprintf('--- regression of mean horizon on log10(dp) ---\n');
fprintf('  levels used                    : %d of %d\n', sum(usable), height(results.levels));
fprintf('  fitted slope s                 : %.4f +- %.4f (1 s.e.)\n', ...
    fit.slope, fit.slopeStandardError);
fprintf('  intercept                      : %.4f\n', fit.intercept);
fprintf('  R^2                            : %.4f\n', fit.rSquared);
fprintf('  lambda = ln10/|s|              : %.4f  [%.4f, %.4f]\n', ...
    results.fit.lambdaFromSlope, results.fit.lambdaLow, results.fit.lambdaHigh);
fprintf('\n  published slope                : %.2f +- %.2f\n', ...
    published.slope, published.slopeUncertainty);
fprintf('  published implies lambda in    : [%.3f, %.3f]\n', ...
    published.lambdaRange(1), published.lambdaRange(2));
fprintf('  <ln|F''|> from the return map    : %.4f\n', returnMapLambda);
fprintf('  which implies slope            : %.4f\n', results.returnMapImpliedSlope);

fprintf('\n--- plane versus accepted crossings ---\n');
fprintf('  plane crossings                : %d\n', results.crossingRatio.planeCrossings);
fprintf('  accepted crossings             : %d\n', results.crossingRatio.acceptedCrossings);
fprintf('  ratio                          : %.4f\n', results.crossingRatio.ratio);
fprintf('  lambda per plane crossing      : %.4f\n', results.lambdaPerPlaneCrossing);
fprintf(['  A map iterate is an accepted crossing, so an exponent measured\n', ...
    '  per iterate exceeds one measured per plane crossing by this ratio.\n']);

fprintf('\n--- classification ---\n');
fprintf('  %s\n', results.classification);
fprintf('  agrees with published fit      : %d\n', agreesWithPublished);
fprintf('  agrees with return-map slope   : %d\n', agreesWithReturnMap);
if results.capBound
    fprintf(['  WARNING: the horizon cap bound in %d trials. Raise\n', ...
        '  options.horizonCap; a binding cap truncates the trend.\n'], ...
        sum(results.levels.trialsAtHorizonCap));
end
if results.anyLevelIncomplete
    fprintf(['  NOTE: not every source succeeded at every dp. The mean is\n', ...
        '  taken over successful sources only, which censors the upper tail\n', ...
        '  and biases the fitted slope toward zero.\n']);
end
fprintf('  runtime                        : %.1f s\n', results.runtimeSeconds);
fprintf('=========================================================\n');
end

% =========================================================================

function fit = local_linear_fit(x, y)
%LOCAL_LINEAR_FIT Least squares with a slope standard error, no toolbox.
x = x(:); y = y(:);
n = numel(x);
X = [ones(n, 1), x];
beta = X \ y;
fit.intercept = beta(1);
fit.slope = beta(2);
residuals = y - X * beta;
fit.residuals = residuals;
if n > 2
    varianceEstimate = sum(residuals.^2) / (n - 2);
    covariance = varianceEstimate * inv(X.' * X); %#ok<MINV>
    fit.slopeStandardError = sqrt(covariance(2, 2));
else
    fit.slopeStandardError = NaN;
end
totalSumSquares = sum((y - mean(y)).^2);
if totalSumSquares > 0
    fit.rSquared = 1 - sum(residuals.^2) / totalSumSquares;
else
    fit.rSquared = NaN;
end
fit.n = n;
end

function ratio = local_crossing_ratio(mapData, params, sourceIndex, horizon)
%LOCAL_CROSSING_RATIO Plane crossings per accepted crossing, uncontrolled.
sourceState = mapData.states(:, sourceIndex);
signatureResult = evaluate_shinbrot_crossing_signature( ...
    sourceState, 0, horizon, params);
ratio.planeCrossings = signatureResult.numberOfPlaneCrossings;
ratio.acceptedCrossings = signatureResult.numberOfAcceptedCrossings;
if ratio.acceptedCrossings > 0
    ratio.ratio = ratio.planeCrossings / ratio.acceptedCrossings;
else
    ratio.ratio = NaN;
end
ratio.horizon = horizon;
ratio.sourceIndex = sourceIndex;
end

function value = local_opt(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end

function config = control_bound_configuration(mapData, params, target, ...
    control, deltaPValues, sourceIndices, randomSeed, horizonCap, ...
    candidateRange, sourceSeparation, returnMapLambda)
config.verifiedBisectionIdentifier = verified_bisection_identifier();
config.bisectionImplementation = verified_bisection_implementation();
config.deterministicDependencyFingerprint = ...
    shinbrot_dependency_fingerprint('deterministic');
config.deltaPValues = deltaPValues(:).';
config.sourceIndices = sourceIndices(:).';
config.randomSeed = randomSeed;
config.horizonCap = horizonCap;
config.candidateRange = candidateRange(:).';
config.sourceSeparation = sourceSeparation;
config.targetX = target.x;
config.targetTolerance = target.tolerance;
config.targetState = target.state(:);
config.sigma = params.sigma;
config.rho = params.rho;
config.beta = params.beta;
config.zSection = params.zSection;
config.xSectionMin = params.xSectionMin;
config.crossingDirection = params.crossingDirection;
config.baseDeltaP = control.deltaP;
config.maxSearchCrossings = horizonCap;
config.bisectionIterations = control.bisectionIterations;
config.maxDiscontinuityIsolationDepth = ...
    get_field_with_default(control, 'maxDiscontinuityIsolationDepth', NaN);
config.maxDiscontinuityIntervals = ...
    get_field_with_default(control, 'maxDiscontinuityIntervals', NaN);
config.returnMapLambda = returnMapLambda;
config.mapDataFingerprint = map_data_fingerprint(mapData);
config.signatureText = stable_config_text(config);
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function fingerprint = map_data_fingerprint(mapData)
fields = {'x', 'eventTimes'};
parts = cell(numel(fields), 1);
for i = 1:numel(fields)
    if isfield(mapData, fields{i})
        values = mapData.(fields{i});
        parts{i} = sprintf('%s:%s', fields{i}, ...
            local_digest(mat2str(values(:).', 17), 'SHA-256'));
    else
        parts{i} = sprintf('%s:missing', fields{i});
    end
end
fingerprint = local_digest(strjoin(parts, newline), 'SHA-256');
end

function textValue = stable_config_text(config)
fields = sort(fieldnames(config));
parts = cell(numel(fields), 1);
for i = 1:numel(fields)
    name = fields{i};
    if strcmp(name, 'signatureText')
        continue;
    end
    value = config.(name);
    if isnumeric(value) || islogical(value)
        valueText = mat2str(value, 17);
    elseif ischar(value)
        valueText = value;
    elseif isstring(value)
        valueText = char(value);
    elseif isstruct(value) && isfield(value, 'aggregateSHA256')
        valueText = value.aggregateSHA256;
    else
        valueText = evalc('disp(value)');
    end
    parts{i} = [name, '=', valueText];
end
textValue = strjoin(parts(~cellfun('isempty', parts)), newline);
end

function digest = local_digest(text, algorithm)
bytes = uint8(text);
try
    engine = java.security.MessageDigest.getInstance(algorithm);
    engine.update(typecast(bytes, 'int8'));
    raw = typecast(engine.digest(), 'uint8');
    digest = lower(reshape(dec2hex(raw, 2).', 1, []));
catch
    weights = mod(1:numel(bytes), 251) + 1;
    digest = sprintf('nojvm%012.0f', ...
        mod(sum(double(bytes) .* weights), 1e12));
end
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
