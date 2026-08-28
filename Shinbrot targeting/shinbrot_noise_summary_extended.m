function extended = shinbrot_noise_summary_extended(trialsTable, calibration, options)
%SHINBROT_NOISE_SUMMARY_EXTENDED Table 4 with intervals and failure mode.
%
%   extended = SHINBROT_NOISE_SUMMARY_EXTENDED(trialsTable, calibration)
%
%   Table 4 currently reports point estimates from ten trials per level and
%   uses none of the diagnostic columns the sweep already records.  Two
%   consequences follow.
%
%   First, a success fraction of 0.2 and one of 0.3 differ by a single trial,
%   yet the table presents them as distinct numbers.  Section 4.1.3 warns
%   against over-reading them in prose; a Wilson score interval makes the
%   warning quantitative and lets the reader see the overlap directly.  The
%   Wilson interval is used rather than the normal approximation because it
%   behaves correctly at zero and at small n, and the sweep contains a level
%   with zero successes.
%
%   Second, the sweep records NumRetargetingFailures but Table 4 reports only
%   NumRetargetingEvents.  Those answer different questions.  A retargeting
%   event that produced no admissible perturbation is a search failure; an
%   event that produced one which the noise then defeated is a control
%   failure.  Without the split, the 80 percent overall failure rate has no
%   attributed mechanism.
%
%   Where a calibration struct from shinbrot_noise_calibration is supplied,
%   the normalised abscissa delta and the ratio sigma/sigma_* are attached to
%   each row, so the table can be read against the published curve.
%
%   trialsTable : the per-trial table, or a path to
%                 noise_trials_bisection_section_retarget.csv
%   calibration : optional output of shinbrot_noise_calibration
%   options     : optional struct
%                   .confidenceLevel  default 0.95
%                   .verbose          default true

if nargin < 2
    calibration = [];
end
if nargin < 3
    options = struct();
end
confidenceLevel = local_option(options, 'confidenceLevel', 0.95);
verbose = local_option(options, 'verbose', true);

if ischar(trialsTable) || isstring(trialsTable)
    trialsTable = readtable(char(trialsTable));
end

sigmaValues = unique(trialsTable.SigmaNoise);
sigmaValues = sort(sigmaValues(:));
z = local_normal_quantile(0.5 * (1 + confidenceLevel));

rows = {};
for i = 1:numel(sigmaValues)
    sigma = sigmaValues(i);
    mask = trialsTable.SigmaNoise == sigma;
    sub = trialsTable(mask, :);

    n = height(sub);
    hits = local_logical(sub.Hit);
    successes = sum(hits);
    fraction = successes / n;
    [ciLow, ciHigh] = local_wilson(successes, n, z);

    successfulCrossings = sub.NumCrossings(hits);
    if isempty(successfulCrossings)
        meanCrossings = NaN;
        medianCrossings = NaN;
    else
        meanCrossings = mean(successfulCrossings);
        medianCrossings = median(successfulCrossings);
    end

    retargetEvents = local_numeric(sub, 'NumRetargetingEvents');
    retargetFailures = local_numeric(sub, 'NumRetargetingFailures');
    rk4Steps = local_numeric(sub, 'NumRk4Steps');
    physicalTime = local_numeric(sub, 'PhysicalTime');

    meanEvents = mean(retargetEvents, 'omitnan');
    meanFailures = mean(retargetFailures, 'omitnan');
    if meanEvents > 0
        failureShare = meanFailures / meanEvents;
    else
        failureShare = NaN;
    end

    rows(end + 1, :) = {sigma, n, successes, fraction, ciLow, ciHigh, ...
        meanCrossings, medianCrossings, ...
        mean(local_numeric(sub, 'FinalTargetError'), 'omitnan'), ...
        std(local_numeric(sub, 'FinalTargetError'), 'omitnan'), ...
        meanEvents, meanFailures, failureShare, ...
        mean(rk4Steps, 'omitnan'), mean(physicalTime, 'omitnan'), ...
        sum(local_logical(local_numeric(sub, 'IntegrationFailed')))}; %#ok<AGROW>
end

extended = cell2table(rows, 'VariableNames', {'sigmaNoise','trials', ...
    'successes','successFraction','wilsonLow','wilsonHigh', ...
    'meanSuccessfulCrossings','medianSuccessfulCrossings', ...
    'meanFinalTargetError','stdFinalTargetError', ...
    'meanRetargetingEvents','meanRetargetingFailures', ...
    'retargetingFailureShare','meanRk4Steps','meanPhysicalTime', ...
    'integrationFailures'});

if ~isempty(calibration) && isstruct(calibration)
    sqrtN = sqrt(calibration.noiseDrawsPerCrossing);
    extended.delta = calibration.kTReferenceDelta * extended.sigmaNoise * sqrtN;
    extended.sigmaOverSigmaStar = extended.sigmaNoise / calibration.sigmaStar;
    extended.aboveThreshold = extended.sigmaNoise > calibration.sigmaStar;
end

if ~verbose
    return;
end

fprintf('\n=== Extended noise summary (%.0f%% Wilson intervals) ===\n', ...
    100 * confidenceLevel);
disp(extended);

nonZero = extended(extended.sigmaNoise > 0, :);
if ~isempty(nonZero)
    fprintf('non-zero-noise trials              : %d\n', sum(nonZero.trials));
    fprintf('non-zero-noise successes           : %d\n', sum(nonZero.successes));
    [lo, hi] = local_wilson(sum(nonZero.successes), sum(nonZero.trials), z);
    fprintf('pooled non-zero success fraction   : %.3f  [%.3f, %.3f]\n', ...
        sum(nonZero.successes) / sum(nonZero.trials), lo, hi);
    overlapping = 0;
    for i = 1:height(nonZero)
        for j = (i + 1):height(nonZero)
            if nonZero.wilsonLow(i) <= nonZero.wilsonHigh(j) && ...
                    nonZero.wilsonLow(j) <= nonZero.wilsonHigh(i)
                overlapping = overlapping + 1;
            end
        end
    end
    totalPairs = height(nonZero) * (height(nonZero) - 1) / 2;
    fprintf('level pairs with overlapping CIs   : %d of %d\n', ...
        overlapping, totalPairs);
    fprintf(['  Overlapping intervals are the quantitative form of the\n', ...
        '  caution already stated in Section 4.1.3: differences between\n', ...
        '  neighbouring success fractions are not resolved by ten trials.\n']);
end

if any(~isnan(extended.meanRetargetingFailures))
    fprintf('\nfailure mechanism\n');
    fprintf(['  retargetingFailureShare is the fraction of retargeting\n', ...
        '  events that returned no admissible perturbation.  A share near\n', ...
        '  zero attributes the trial failures to the noise defeating a\n', ...
        '  valid control; a large share attributes them to the search\n', ...
        '  itself finding nothing.\n']);
end
fprintf('=========================================================\n');
end

% -------------------------------------------------------------------------

function [lo, hi] = local_wilson(successes, n, z)
%LOCAL_WILSON Wilson score interval, valid at zero successes and small n.
if n == 0
    lo = NaN;
    hi = NaN;
    return;
end
p = successes / n;
denominator = 1 + z^2 / n;
centre = (p + z^2 / (2 * n)) / denominator;
halfWidth = (z / denominator) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2));
lo = max(0, centre - halfWidth);
hi = min(1, centre + halfWidth);
end

function q = local_normal_quantile(p)
%LOCAL_NORMAL_QUANTILE Standard normal inverse CDF without the Statistics Toolbox.
q = -sqrt(2) * erfcinv(2 * p);
end

function v = local_numeric(tbl, name)
if ismember(name, tbl.Properties.VariableNames)
    v = tbl.(name);
    if iscell(v)
        v = cellfun(@local_to_double, v);
    end
    v = double(v);
else
    v = NaN(height(tbl), 1);
end
end

function d = local_to_double(value)
if ischar(value) || isstring(value)
    d = str2double(value);
else
    d = double(value);
end
end

function tf = local_logical(v)
if islogical(v)
    tf = v;
else
    tf = double(v) ~= 0;
end
tf = tf(:);
end

function value = local_option(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
