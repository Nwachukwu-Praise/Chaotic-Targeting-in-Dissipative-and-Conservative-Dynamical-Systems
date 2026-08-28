function calibration = shinbrot_noise_calibration(mapData, target, dt, ...
    sigmaNoiseValues, options)
%SHINBROT_NOISE_CALIBRATION Place the noise sweep on the published abscissa.
%
%   calibration = SHINBROT_NOISE_CALIBRATION(mapData, target, dt, sigmaNoiseValues)
%
%   Subsection 3.1.9 promises three quantities in Section 4.1.3 that the
%   draft does not currently report:
%
%     1. Tbar, the mean interval between accepted section crossings of the
%        nominal map, and N = Tbar/dt, the number of noise draws applied per
%        accepted crossing;
%     2. the conversion delta ~ 0.572 * sigma * sqrt(N) of equation (23),
%        which maps the raw per-step, per-coordinate standard deviation used
%        here onto the normalised abscissa of the published noise-response
%        curve;
%     3. sigma_* = eps_t / sqrt(N) of equation (24), the level above which
%        the displacement accumulated between a control decision and the
%        next crossing exceeds the target half-width.
%
%   Without (3) the reader cannot tell whether a low success fraction is a
%   property of the method or a property of the levels that were tested.
%
%   Inputs
%     mapData           : output of generate_return_map for the NOMINAL map
%                         (pControl = 0).  Only mapData.eventTimes is used.
%     target            : struct with field tolerance (eps_t).
%     dt                : fixed integration step used in the noise study.
%     sigmaNoiseValues  : the levels of equation (22).
%     options           : optional struct.
%       .burnIn                  crossings discarded as transient (default 100)
%       .kTReferenceDelta        plotted level at which kT = 1 (default 0.572)
%       .adaptiveFinalError      Table 2 final target error, adaptive solver
%       .fixedStepFinalError     Table 4 sigma = 0 final target error
%       .verbose                 print the report (default true)
%
%   Nothing here re-runs the sweep.  It is a reinterpretation of results you
%   already have, plus one measurement on the nominal map.

if nargin < 5
    options = struct();
end
burnIn = local_option(options, 'burnIn', 100);
kTReferenceDelta = local_option(options, 'kTReferenceDelta', 0.572);
verbose = local_option(options, 'verbose', true);

if ~isfield(mapData, 'eventTimes') || isempty(mapData.eventTimes)
    error('shinbrot_noise_calibration:NoEventTimes', ...
        ['mapData.eventTimes is required.  Regenerate the nominal map with ', ...
        'generate_return_map, which stores it.']);
end

% ---------------------------------------------------------------- Tbar
eventTimes = mapData.eventTimes(:).';
eventTimes = eventTimes(isfinite(eventTimes));
if numel(eventTimes) < burnIn + 10
    burnIn = 0;
end
returnTimes = diff(eventTimes(burnIn + 1:end));
returnTimes = returnTimes(isfinite(returnTimes) & returnTimes > 0);
if isempty(returnTimes)
    error('shinbrot_noise_calibration:NoReturnTimes', ...
        'No usable inter-crossing intervals were found.');
end

Tbar = mean(returnTimes);
N = Tbar / dt;
sqrtN = sqrt(N);
epsilonT = target.tolerance;
sigmaStar = epsilonT / sqrtN;

calibration.burnInCrossings = burnIn;
calibration.numberOfReturnIntervals = numel(returnTimes);
calibration.meanReturnTime = Tbar;
calibration.medianReturnTime = median(returnTimes);
calibration.stdReturnTime = std(returnTimes);
calibration.minReturnTime = min(returnTimes);
calibration.maxReturnTime = max(returnTimes);
calibration.integrationStep = dt;
calibration.noiseDrawsPerCrossing = N;
calibration.targetTolerance = epsilonT;
calibration.sigmaStar = sigmaStar;
calibration.kTReferenceDelta = kTReferenceDelta;
% delta = 0.572 * sigma * sqrt(N); kT = 1 is plotted at delta = 0.572, so
% the sigma matching the published unit-noise reference is 1/sqrt(N).
calibration.sigmaAtUnitKT = 1 / sqrtN;

% ------------------------------------------------------- per-level table
sigmaNoiseValues = sigmaNoiseValues(:);
delta = kTReferenceDelta * sigmaNoiseValues * sqrtN;
ratioToThreshold = sigmaNoiseValues / sigmaStar;
perCrossingDisplacement = sigmaNoiseValues * sqrtN;
belowThreshold = sigmaNoiseValues <= sigmaStar;
withinPublishedAbscissa = delta <= 1;

calibration.levels = table(sigmaNoiseValues, delta, ratioToThreshold, ...
    perCrossingDisplacement, belowThreshold, withinPublishedAbscissa, ...
    'VariableNames', {'sigmaNoise','delta','sigmaOverSigmaStar', ...
    'perCrossingDisplacement','belowThreshold','deltaAtMostOne'});

% --------------------------------------- adaptive versus fixed-step check
adaptiveError = local_option(options, 'adaptiveFinalError', NaN);
fixedStepError = local_option(options, 'fixedStepFinalError', NaN);
calibration.adaptiveFinalError = adaptiveError;
calibration.fixedStepFinalError = fixedStepError;
calibration.discretisationDifference = abs(adaptiveError - fixedStepError);
calibration.discretisationRatio = adaptiveError / fixedStepError;
calibration.discretisationAsFractionOfTolerance = ...
    calibration.discretisationDifference / epsilonT;

if ~verbose
    return;
end

fprintf('\n=== Shinbrot noise calibration (equations 23 and 24) ===\n');
fprintf('nominal map, transient discarded : %d crossings\n', burnIn);
fprintf('inter-crossing intervals used    : %d\n', numel(returnTimes));
fprintf('mean return time  Tbar           : %.6g\n', Tbar);
fprintf('  median / std                   : %.6g / %.6g\n', ...
    calibration.medianReturnTime, calibration.stdReturnTime);
fprintf('  min / max                      : %.6g / %.6g\n', ...
    calibration.minReturnTime, calibration.maxReturnTime);
fprintf('integration step  dt             : %.6g\n', dt);
fprintf('noise draws per crossing  N      : %.6g   (sqrt(N) = %.4g)\n', N, sqrtN);
fprintf('target half-width eps_t          : %.6g\n', epsilonT);
fprintf('threshold  sigma_* = eps_t/sqrtN : %.4g\n', sigmaStar);
fprintf('sigma giving the published kT=1  : %.4g   (delta = %.3g)\n', ...
    calibration.sigmaAtUnitKT, kTReferenceDelta);
fprintf('\nper-level conversion\n');
disp(calibration.levels);

nonZero = sigmaNoiseValues(sigmaNoiseValues > 0);
if ~isempty(nonZero)
    fprintf('smallest non-zero level tested   : %.4g\n', min(nonZero));
    fprintf('  as a multiple of sigma_*       : %.4gx\n', min(nonZero) / sigmaStar);
    fprintf('levels at or below sigma_*       : %d of %d\n', ...
        sum(belowThreshold & sigmaNoiseValues > 0), numel(nonZero));
    fprintf('levels with delta <= 1           : %d of %d\n', ...
        sum(withinPublishedAbscissa & sigmaNoiseValues > 0), numel(nonZero));
end

if isfinite(calibration.discretisationDifference)
    fprintf('\nadaptive versus fixed-step discretisation (subsection 3.1.9)\n');
    fprintf('  adaptive deterministic error   : %.6g\n', adaptiveError);
    fprintf('  fixed-step sigma = 0 error     : %.6g\n', fixedStepError);
    fprintf('  absolute difference            : %.6g\n', ...
        calibration.discretisationDifference);
    fprintf('  ratio                          : %.4g\n', calibration.discretisationRatio);
    fprintf('  as a fraction of eps_t         : %.4g\n', ...
        calibration.discretisationAsFractionOfTolerance);
else
    fprintf(['\nadaptive/fixed-step comparison skipped: supply ', ...
        'options.adaptiveFinalError and options.fixedStepFinalError.\n']);
end

fprintf(['\nReading: the sweep samples the descending branch of the noise ', ...
    'response\nwherever sigma exceeds sigma_*, and lies off the published ', ...
    'abscissa\nwherever delta exceeds one.  Neither is a property of the ', ...
    'targeting method.\n']);
fprintf('=========================================================\n');
end

function value = local_option(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
