function diagnostics = shinbrot_return_map_diagnostics(mapData, target, options)
%SHINBROT_RETURN_MAP_DIAGNOSTICS Slope exponent and one-dimensionality.
%
%   diagnostics = SHINBROT_RETURN_MAP_DIAGNOSTICS(mapData, target)
%
%   Two properties of the nominal return map that the draft asserts but does
%   not measure.
%
%   (a) The mean logarithmic slope of equation (27),
%
%           <ln|F'|> = (1/M) sum_j ln| (X_{n+1}^{(j)} - X_{n+1}^{(i)}) /
%                                      (X_n^{(j)}   - X_n^{(i)}  ) |,
%
%       estimated over pairs of crossings adjacent in X_n and lying on the
%       same branch.  Subsections 2.1.1 and 3.1.10 note that Shinbrot's
%       lambda is not uniquely fixed and offer this as a third candidate
%       alongside the largest Lyapunov exponent and the topological entropy.
%       The implied regression slope s = -ln(10)/lambda of equation (26) is
%       reported so it can be compared with the published -2.49 +- 0.09.
%
%   (b) The thickness of the return map.  Subsection 3.1.2 records that F is
%       single-valued in X_n only to the accuracy of the one-dimensional
%       approximation, the attractor having fractal dimension slightly above
%       two.  That is the modelling assumption on which the entire scalar
%       formulation rests, and it is currently only asserted.  Here it is
%       quantified: within sliding windows in X_n a local straight line is
%       fitted and the residual spread in X_{n+1} is measured.  The number
%       that matters is the ratio of that spread to the target half-width
%       eps_t, because the method only needs one-dimensionality to hold at
%       the scale of the target.
%
%   options
%     .burnIn            crossings discarded as transient (default 100)
%     .branchJumpFactor  a pair is treated as crossing a branch boundary when
%                        |dX_{n+1}| exceeds this multiple of the median
%                        |dX_{n+1}| over adjacent pairs (default 8)
%     .minSeparation     pairs closer than this in X_n are discarded as
%                        numerically unreliable (default 1e-6)
%     .windowCount       number of sliding windows for the thickness
%                        estimate (default 40)
%     .minWindowPoints   minimum crossings per window (default 12)
%     .verbose           print the report (default true)

if nargin < 3
    options = struct();
end
burnIn = local_option(options, 'burnIn', 100);
branchJumpFactor = local_option(options, 'branchJumpFactor', 8);
minSeparation = local_option(options, 'minSeparation', 1e-6);
windowCount = local_option(options, 'windowCount', 40);
minWindowPoints = local_option(options, 'minWindowPoints', 12);
verbose = local_option(options, 'verbose', true);

x = mapData.x(:).';
if numel(x) < burnIn + 50
    burnIn = 0;
end
x = x(burnIn + 1:end);
Xn = x(1:end - 1);
Xnext = x(2:end);

diagnostics.burnInCrossings = burnIn;
diagnostics.numberOfPairs = numel(Xn);
diagnostics.targetTolerance = target.tolerance;

% ------------------------------------------------- (a) mean log slope
[sorted, order] = sort(Xn);
sortedNext = Xnext(order);

dX = diff(sorted);
dXnext = diff(sortedNext);

usable = dX > minSeparation;
medianStep = median(abs(dXnext(usable)));
sameBranch = usable & abs(dXnext) <= branchJumpFactor * medianStep;

slopes = dXnext(sameBranch) ./ dX(sameBranch);
logSlopes = log(abs(slopes));
logSlopes = logSlopes(isfinite(logSlopes));

diagnostics.branchJumpFactor = branchJumpFactor;
diagnostics.medianAdjacentImageStep = medianStep;
diagnostics.adjacentPairsConsidered = sum(usable);
diagnostics.adjacentPairsSameBranch = sum(sameBranch);
diagnostics.adjacentPairsRejectedAsBranchJump = sum(usable & ~sameBranch);
diagnostics.meanLogSlope = mean(logSlopes);
diagnostics.medianLogSlope = median(logSlopes);
diagnostics.stdLogSlope = std(logSlopes);
diagnostics.impliedRegressionSlope = -log(10) / diagnostics.meanLogSlope;
diagnostics.publishedRegressionSlope = -2.49;
diagnostics.publishedRegressionSlopeUncertainty = 0.09;
diagnostics.publishedLambdaRange = ...
    log(10) ./ [2.49 + 0.09, 2.49 - 0.09];

% Sensitivity: the same estimate without any branch rejection, so the
% dependence on the declared criterion is visible rather than hidden.
allSlopes = log(abs(dXnext(usable) ./ dX(usable)));
allSlopes = allSlopes(isfinite(allSlopes));
diagnostics.meanLogSlopeNoBranchRejection = mean(allSlopes);

% ------------------------------------------------- (b) map thickness
lo = min(Xn);
hi = max(Xn);
edges = linspace(lo, hi, windowCount + 1);
rows = {};
allResiduals = [];
for w = 1:windowCount
    inWindow = Xn >= edges(w) & Xn < edges(w + 1);
    if w == windowCount
        inWindow = inWindow | Xn == edges(end);
    end
    if sum(inWindow) < minWindowPoints
        continue;
    end
    wx = Xn(inWindow).';
    wy = Xnext(inWindow).';
    % Reject windows that straddle a branch boundary: a straight line is not
    % a meaningful local model there, and the residual would measure the
    % branch gap rather than the map thickness.
    spread = max(wy) - min(wy);
    localMedianStep = median(abs(diff(sort(wy))));
    if spread > branchJumpFactor * localMedianStep * numel(wy)
        continue;
    end
    A = [ones(numel(wx), 1), wx];
    coefficients = A \ wy;
    residuals = wy - A * coefficients;
    allResiduals = [allResiduals; residuals]; %#ok<AGROW>
    rows(end + 1, :) = {w, 0.5 * (edges(w) + edges(w + 1)), numel(wx), ...
        coefficients(2), sqrt(mean(residuals.^2)), max(abs(residuals))}; %#ok<AGROW>
end

if isempty(rows)
    diagnostics.windows = table();
    diagnostics.thicknessRms = NaN;
    diagnostics.thicknessMax = NaN;
else
    diagnostics.windows = cell2table(rows, 'VariableNames', ...
        {'window','centreXn','points','localSlope','rmsResidual','maxResidual'});
    diagnostics.thicknessRms = sqrt(mean(allResiduals.^2));
    diagnostics.thicknessMax = max(abs(allResiduals));
end
diagnostics.thicknessRmsOverTolerance = ...
    diagnostics.thicknessRms / target.tolerance;
diagnostics.thicknessMaxOverTolerance = ...
    diagnostics.thicknessMax / target.tolerance;

if ~verbose
    return;
end

fprintf('\n=== Nominal return-map diagnostics ===\n');
fprintf('crossings used (transient %d discarded) : %d pairs\n', ...
    burnIn, diagnostics.numberOfPairs);

fprintf('\n(a) mean logarithmic slope, equation (27)\n');
fprintf('  adjacent pairs considered            : %d\n', ...
    diagnostics.adjacentPairsConsidered);
fprintf('  retained as same-branch              : %d\n', ...
    diagnostics.adjacentPairsSameBranch);
fprintf('  rejected as branch jumps             : %d\n', ...
    diagnostics.adjacentPairsRejectedAsBranchJump);
fprintf('  <ln|F''|>                             : %.4f\n', diagnostics.meanLogSlope);
fprintf('  median / std                         : %.4f / %.4f\n', ...
    diagnostics.medianLogSlope, diagnostics.stdLogSlope);
fprintf('  without branch rejection             : %.4f  (sensitivity check)\n', ...
    diagnostics.meanLogSlopeNoBranchRejection);
fprintf('  implied slope s = -ln10/lambda       : %.4f\n', ...
    diagnostics.impliedRegressionSlope);
fprintf('  published s                          : %.2f +- %.2f\n', ...
    diagnostics.publishedRegressionSlope, ...
    diagnostics.publishedRegressionSlopeUncertainty);
fprintf('  published s implies lambda in        : [%.3f, %.3f]\n', ...
    diagnostics.publishedLambdaRange(1), diagnostics.publishedLambdaRange(2));

fprintf('\n(b) one-dimensionality of the return map, subsection 3.1.2\n');
fprintf('  windows accepted                     : %d of %d\n', ...
    height(diagnostics.windows), windowCount);
fprintf('  RMS residual in X_{n+1}              : %.4g\n', diagnostics.thicknessRms);
fprintf('  max residual                         : %.4g\n', diagnostics.thicknessMax);
fprintf('  target half-width eps_t              : %.4g\n', target.tolerance);
fprintf('  RMS residual / eps_t                 : %.4g\n', ...
    diagnostics.thicknessRmsOverTolerance);
fprintf('  max residual / eps_t                 : %.4g\n', ...
    diagnostics.thicknessMaxOverTolerance);
if diagnostics.thicknessRmsOverTolerance < 1
    fprintf(['  The scalar map is single-valued to better than the target\n', ...
        '  half-width, so the one-dimensional reduction is adequate at\n', ...
        '  the scale the targeting problem actually uses.\n']);
else
    fprintf(['  The residual spread is comparable to or larger than the\n', ...
        '  target half-width.  The scalar reduction is then a real source\n', ...
        '  of error in the targeting calculation and should be reported\n', ...
        '  as such rather than assumed away.\n']);
end
fprintf('======================================\n');
end

function value = local_option(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
