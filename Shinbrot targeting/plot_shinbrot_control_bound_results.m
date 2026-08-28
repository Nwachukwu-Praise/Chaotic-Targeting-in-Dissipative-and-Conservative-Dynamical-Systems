function plot_shinbrot_control_bound_results(results)
%PLOT_SHINBROT_CONTROL_BOUND_RESULTS Plot repaired paper-bisection scaling.

validate_control_bound_result(results);

levels = results.levels;
usable = levels.successes > 0 & isfinite(levels.meanHorizon);
x = levels.log10DeltaP(usable);
y = levels.meanHorizon(usable);
se = levels.stdHorizon(usable) ./ sqrt(max(levels.successes(usable), 1));
anchorX = mean(x);
anchorY = polyval([results.fit.slope, results.fit.intercept], anchorX);
xLine = linspace(min(x), max(x), 100);
fitLine = polyval([results.fit.slope, results.fit.intercept], xLine);
publishedLine = anchorY + results.published.slope * (xLine - anchorX);
returnMapLine = anchorY + results.returnMapImpliedSlope * (xLine - anchorX);

figure('Name', 'Scaling in the control bound', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
errorbar(x, y, se, 'o', 'LineWidth', 1.2, ...
    'MarkerFaceColor', [0.10, 0.35, 0.65], ...
    'Color', [0.10, 0.35, 0.65], ...
    'DisplayName', 'mean horizon +- s.e.');
hold on;
plot(xLine, fitLine, 'k-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('fit slope %.3g', results.fit.slope));
plot(xLine, publishedLine, '--', 'LineWidth', 1.3, ...
    'Color', [0.75, 0.25, 0.10], ...
    'DisplayName', 'published slope');
plot(xLine, returnMapLine, ':', 'LineWidth', 1.6, ...
    'Color', [0.10, 0.55, 0.30], ...
    'DisplayName', 'return-map slope');
grid on;
xlabel('log_{10}(\Delta p)');
ylabel('mean targeting horizon n_t');
title('Mean horizon versus control bound');
legend('Location', 'best');
hold off;

nexttile;
fitLambda = results.fit.lambdaFromSlope;
lambdaLow = min(results.fit.lambdaLow, results.fit.lambdaHigh);
lambdaHigh = max(results.fit.lambdaLow, results.fit.lambdaHigh);
errorbar(1, fitLambda, fitLambda - lambdaLow, lambdaHigh - fitLambda, ...
    'o', 'LineWidth', 1.5, 'MarkerFaceColor', [0.10, 0.35, 0.65], ...
    'Color', [0.10, 0.35, 0.65], 'DisplayName', 'fit');
hold on;
yline(results.returnMapLambda, ':', 'Color', [0.10, 0.55, 0.30], ...
    'LineWidth', 1.6, 'DisplayName', 'return-map \lambda');
publishedLow = min(results.published.lambdaRange);
publishedHigh = max(results.published.lambdaRange);
patch([0.65, 1.35, 1.35, 0.65], ...
    [publishedLow, publishedLow, publishedHigh, publishedHigh], ...
    [0.75, 0.25, 0.10], 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'DisplayName', 'published range');
grid on;
xlim([0.5, 1.5]);
set(gca, 'XTick', 1, 'XTickLabel', {'\lambda estimate'});
ylabel('\lambda');
title('Exponent comparison');
legend('Location', 'best');
hold off;

fprintf(['Control-bound plot uses repaired %s results with %d deltaP ', ...
    'levels and %d source states.\n'], ...
    results.verifiedBisectionIdentifier, height(levels), ...
    numel(results.sourceIndices));
end

function validate_control_bound_result(results)
required = {'levels', 'fit', 'published', 'returnMapLambda', ...
    'returnMapImpliedSlope', 'sourceIndices', ...
    'verifiedBisectionIdentifier', 'configuration', 'trials'};
if ~isstruct(results) || ~all(isfield(results, required))
    error('Control-bound result lacks required repaired-result fields.');
end
verifiedIdentifier = verified_bisection_identifier();
if ~strcmp(results.verifiedBisectionIdentifier, verifiedIdentifier)
    error('Control-bound result does not use the verified paper bisection identifier.');
end
if ~isfield(results.configuration, 'signatureText')
    error('Control-bound result lacks a configuration signature.');
end
if ~istable(results.trials) || ...
        ~any(strcmp(results.trials.Properties.VariableNames, 'searchMethod')) || ...
        ~all(strcmp(results.trials.searchMethod, verifiedIdentifier))
    error('Control-bound trials do not all report the verified bisection method.');
end
end
