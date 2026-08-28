function figures = plot_shinbrot_source_ensemble( ...
    pairedResultsTable, outputDirectory)
%PLOT_SHINBROT_SOURCE_ENSEMBLE Publication-quality paired comparisons.
%
% figures = plot_shinbrot_source_ensemble(pairedResultsTable)
% figures = plot_shinbrot_source_ensemble(pairedResultsTable, outputDirectory)
%
% The observations are discrete source states. Consequently, this function
% does not join one source state to the next. Light connectors are used only
% to pair the two methods evaluated from the same source state.
%
% If outputDirectory is supplied, the four figures are also exported as
% vector SVG files. The one-input call remains fully compatible with the
% original function.

if nargin < 2
    outputDirectory = '';
end
if isstring(outputDirectory) && isscalar(outputDirectory)
    outputDirectory = char(outputDirectory);
end

required = {'SourceIndex', 'ShinbrotReplayHit', 'HybridReplayHit', ...
    'ShinbrotAbsoluteP', 'HybridAbsoluteP', 'AbsolutePDifference', ...
    'ShinbrotHorizon', 'HybridHorizon', 'SameHorizon', 'BothReplayHit'};
if ~all(ismember(required, pairedResultsTable.Properties.VariableNames))
    error('pairedResultsTable does not contain the required variables.');
end
if isempty(pairedResultsTable)
    error('pairedResultsTable is empty.');
end

% Plot the source states in ascending return-map order.
[sourceIndex, order] = sort(pairedResultsTable.SourceIndex(:));
paired = pairedResultsTable(order, :);

shinbrotVerified = logical(paired.ShinbrotReplayHit);
hybridVerified = logical(paired.HybridReplayHit);
bothVerified = logical(paired.BothReplayHit);

shinbrotMagnitude = paired.ShinbrotAbsoluteP;
hybridMagnitude = paired.HybridAbsoluteP;
pairedDifference = paired.AbsolutePDifference;
shinbrotHorizon = paired.ShinbrotHorizon;
hybridHorizon = paired.HybridHorizon;

% Colour-blind-safe, print-friendly palette.
colour.shinbrot = [0.1216, 0.3059, 0.4745];
colour.hybrid = [0.8353, 0.3686, 0.0000];
colour.equal = [0.35, 0.37, 0.40];
colour.connector = [0.76, 0.78, 0.81];
colour.differentHorizon = [0.62, 0.12, 0.18];

numberSources = height(paired);
sourceOffset = source_marker_offset(sourceIndex);
figures = gobjects(4, 1);
axesHandles = gobjects(4, 1);

%% 1. Returned perturbation magnitude: paired dumbbell plot
[figures(1), axesHandles(1)] = new_ensemble_figure( ...
    'Returned perturbation magnitude', [100, 100, 1040, 600]);
ax = axesHandles(1);

pairMask = bothVerified & isfinite(shinbrotMagnitude) & ...
    isfinite(hybridMagnitude);
for row = find(pairMask).'
    plot(ax, sourceIndex(row) + [-sourceOffset, sourceOffset], ...
        [shinbrotMagnitude(row), hybridMagnitude(row)], '-', ...
        'Color', colour.connector, 'LineWidth', 0.85, ...
        'HandleVisibility', 'off');
end

shinbrotMask = shinbrotVerified & isfinite(shinbrotMagnitude);
hybridMask = hybridVerified & isfinite(hybridMagnitude);
plot(ax, sourceIndex(shinbrotMask) - sourceOffset, ...
    shinbrotMagnitude(shinbrotMask), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', colour.shinbrot, 'LineWidth', 1.45, ...
    'DisplayName', 'Shinbrot bisection');
plot(ax, sourceIndex(hybridMask) + sourceOffset, ...
    hybridMagnitude(hybridMask), 'd', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', colour.hybrid, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', 'Hybrid grid-bisection');

plot_failed_magnitudes(ax, paired, 'Shinbrot', ...
    -sourceOffset, colour.shinbrot);
plot_failed_magnitudes(ax, paired, 'Hybrid', ...
    sourceOffset, colour.hybrid);

set_source_limits(ax, sourceIndex, sourceOffset);
allVerifiedMagnitudes = [shinbrotMagnitude(shinbrotMask); ...
    hybridMagnitude(hybridMask)];
set_nonnegative_limits(ax, allVerifiedMagnitudes);
xlabel(ax, 'Source return-map index');
ylabel(ax, '$|p_{\mathrm{control}}|$', 'Interpreter', 'latex');
title(ax, sprintf( ...
    'Returned perturbation magnitude across %d source states', ...
    numberSources));
show_legend(ax, 'northoutside', 'horizontal');

%% 2. Signed within-source perturbation difference
[figures(2), axesHandles(2)] = new_ensemble_figure( ...
    'Paired perturbation difference', [120, 120, 1040, 600]);
ax = axesHandles(2);

comparisonMask = bothVerified & isfinite(pairedDifference);
negativeMask = comparisonMask & pairedDifference < 0;
positiveMask = comparisonMask & pairedDifference > 0;
equalMask = comparisonMask & pairedDifference == 0;

yline(ax, 0, '-', 'Color', colour.equal, 'LineWidth', 1.0, ...
    'HandleVisibility', 'off');
for row = find(comparisonMask).'
    plot(ax, [sourceIndex(row), sourceIndex(row)], ...
        [0, pairedDifference(row)], '-', ...
        'Color', colour.connector, 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
end

plot(ax, sourceIndex(negativeMask), pairedDifference(negativeMask), 'v', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', colour.hybrid, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Hybrid smaller (n = %d)', sum(negativeMask)));
plot(ax, sourceIndex(positiveMask), pairedDifference(positiveMask), '^', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', colour.shinbrot, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Shinbrot smaller (n = %d)', sum(positiveMask)));
plot(ax, sourceIndex(equalMask), pairedDifference(equalMask), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', colour.equal, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Equal magnitude (n = %d)', sum(equalMask)));

set_source_limits(ax, sourceIndex, 0);
set_symmetric_limits(ax, pairedDifference(comparisonMask));
xlabel(ax, 'Source return-map index');
ylabel(ax, ...
    '$\Delta |p|=|p_{\mathrm{hybrid}}|-|p_{\mathrm{Shinbrot}}|$', ...
    'Interpreter', 'latex');
title(ax, {'Within-source perturbation-magnitude difference', ...
    'Negative values favour the hybrid method'});
show_legend(ax, 'northoutside', 'horizontal');

%% 3. Selected targeting horizon: offset paired markers
[figures(3), axesHandles(3)] = new_ensemble_figure( ...
    'Selected targeting horizon', [140, 140, 1040, 600]);
ax = axesHandles(3);

shinbrotHorizonMask = shinbrotVerified & isfinite(shinbrotHorizon);
hybridHorizonMask = hybridVerified & isfinite(hybridHorizon);
horizonPairMask = bothVerified & isfinite(shinbrotHorizon) & ...
    isfinite(hybridHorizon);
differentHorizon = horizonPairMask & ...
    shinbrotHorizon ~= hybridHorizon;

differentLabelUsed = false;
for row = find(horizonPairMask).'
    if differentHorizon(row)
        lineColour = colour.differentHorizon;
        lineWidth = 1.6;
        if ~differentLabelUsed
            displayName = 'Different selected horizon';
            handleVisibility = 'on';
            differentLabelUsed = true;
        else
            displayName = '';
            handleVisibility = 'off';
        end
    else
        lineColour = colour.connector;
        lineWidth = 0.8;
        displayName = '';
        handleVisibility = 'off';
    end
    plot(ax, sourceIndex(row) + [-sourceOffset, sourceOffset], ...
        [shinbrotHorizon(row), hybridHorizon(row)], '-', ...
        'Color', lineColour, 'LineWidth', lineWidth, ...
        'DisplayName', displayName, ...
        'HandleVisibility', handleVisibility);
end

plot(ax, sourceIndex(shinbrotHorizonMask) - sourceOffset, ...
    shinbrotHorizon(shinbrotHorizonMask), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', colour.shinbrot, 'LineWidth', 1.45, ...
    'DisplayName', 'Shinbrot bisection');
plot(ax, sourceIndex(hybridHorizonMask) + sourceOffset, ...
    hybridHorizon(hybridHorizonMask), 'd', ...
    'LineStyle', 'none', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', colour.hybrid, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', 'Hybrid grid-bisection');

set_source_limits(ax, sourceIndex, sourceOffset);
set_horizon_limits(ax, [shinbrotHorizon(shinbrotHorizonMask); ...
    hybridHorizon(hybridHorizonMask)]);
xlabel(ax, 'Source return-map index');
ylabel(ax, 'Selected targeting horizon, $n$', 'Interpreter', 'latex');
title(ax, 'Selected targeting horizon across source states');
show_legend(ax, 'northoutside', 'horizontal');

%% 4. Direct paired comparison against the equality line
[figures(4), axesHandles(4)] = new_ensemble_figure( ...
    'Paired perturbation comparison', [160, 80, 760, 690]);
ax = axesHandles(4);
ax.XGrid = 'on';

validX = shinbrotMagnitude(comparisonMask);
validY = hybridMagnitude(comparisonMask);
allMagnitudes = [validX; validY];
comparisonLimits = nonnegative_comparison_limits(allMagnitudes);
plot(ax, comparisonLimits, comparisonLimits, '--', ...
    'Color', colour.equal, 'LineWidth', 1.1, ...
    'DisplayName', 'Equal magnitude');

plot(ax, shinbrotMagnitude(negativeMask), ...
    hybridMagnitude(negativeMask), 'd', ...
    'LineStyle', 'none', 'MarkerSize', 7.0, ...
    'MarkerFaceColor', colour.hybrid, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Hybrid smaller (n = %d)', sum(negativeMask)));
plot(ax, shinbrotMagnitude(positiveMask), ...
    hybridMagnitude(positiveMask), 'o', ...
    'LineStyle', 'none', 'MarkerSize', 7.0, ...
    'MarkerFaceColor', colour.shinbrot, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Shinbrot smaller (n = %d)', sum(positiveMask)));
plot(ax, shinbrotMagnitude(equalMask), hybridMagnitude(equalMask), 's', ...
    'LineStyle', 'none', 'MarkerSize', 6.0, ...
    'MarkerFaceColor', colour.equal, ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.75, ...
    'DisplayName', sprintf('Equal magnitude (n = %d)', sum(equalMask)));

xlim(ax, comparisonLimits);
ylim(ax, comparisonLimits);
axis(ax, 'square');
xlabel(ax, 'Shinbrot bisection, $|p_{\mathrm{control}}|$', ...
    'Interpreter', 'latex');
ylabel(ax, 'Hybrid grid-bisection, $|p_{\mathrm{control}}|$', ...
    'Interpreter', 'latex');
title(ax, 'Paired perturbation magnitudes');
show_legend(ax, 'northoutside', 'vertical');

if ~isempty(outputDirectory)
    export_ensemble_svgs(figures, outputDirectory);
end
end

function [figureHandle, axesHandle] = new_ensemble_figure(name, position)
%NEW_ENSEMBLE_FIGURE Create one consistently formatted figure and axes.

figureHandle = figure('Name', name, 'NumberTitle', 'off', ...
    'Color', 'white', 'Renderer', 'painters', ...
    'Units', 'pixels', 'Position', position);
axesHandle = axes('Parent', figureHandle);
hold(axesHandle, 'on');

axesHandle.FontName = 'Helvetica';
axesHandle.FontSize = 11;
axesHandle.LineWidth = 0.85;
axesHandle.TickDir = 'out';
axesHandle.TickLength = [0.012, 0.012];
axesHandle.Box = 'off';
axesHandle.Layer = 'top';
axesHandle.Color = 'white';
axesHandle.GridColor = [0.84, 0.86, 0.89];
axesHandle.GridAlpha = 0.70;
axesHandle.XGrid = 'off';
axesHandle.YGrid = 'on';
end

function offset = source_marker_offset(sourceIndex)
%SOURCE_MARKER_OFFSET Separate paired markers without obscuring the index.

sourceRange = max(sourceIndex) - min(sourceIndex);
offset = max(0.25, 0.003 * sourceRange);
end

function set_source_limits(ax, sourceIndex, markerOffset)
%SET_SOURCE_LIMITS Add restrained horizontal padding around source indices.

sourceMinimum = min(sourceIndex);
sourceMaximum = max(sourceIndex);
sourceRange = sourceMaximum - sourceMinimum;
padding = max(1, 0.025 * max(sourceRange, 1)) + markerOffset;
xlim(ax, [sourceMinimum - padding, sourceMaximum + padding]);
end

function set_nonnegative_limits(ax, values)
%SET_NONNEGATIVE_LIMITS Set a zero-based vertical scale with headroom.

values = values(isfinite(values));
if isempty(values)
    ylim(ax, [0, 1]);
    return;
end
maximumValue = max(values);
if maximumValue <= 0
    upperLimit = 1;
else
    upperLimit = maximumValue + max(0.005, 0.05 * maximumValue);
end
ylim(ax, [0, upperLimit]);
end

function set_symmetric_limits(ax, values)
%SET_SYMMETRIC_LIMITS Centre a paired-difference scale on zero.

values = values(isfinite(values));
if isempty(values)
    limit = 1;
else
    limit = max(abs(values));
    if limit <= 0
        limit = 1;
    else
        limit = 1.08 * limit;
    end
end
ylim(ax, [-limit, limit]);
end

function set_horizon_limits(ax, values)
%SET_HORIZON_LIMITS Use integer ticks for the discrete targeting horizon.

values = values(isfinite(values));
if isempty(values)
    ylim(ax, [0, 1]);
    yticks(ax, 0:1);
    return;
end
lowerLimit = max(0, floor(min(values)) - 1);
upperLimit = ceil(max(values)) + 1;
if upperLimit <= lowerLimit
    upperLimit = lowerLimit + 1;
end
ylim(ax, [lowerLimit, upperLimit]);
yticks(ax, lowerLimit:upperLimit);
end

function limits = nonnegative_comparison_limits(values)
%NONNEGATIVE_COMPARISON_LIMITS Return equal limits for a paired scatter plot.

values = values(isfinite(values));
if isempty(values)
    limits = [0, 1];
    return;
end
maximumValue = max(values);
if maximumValue <= 0
    upperLimit = 1;
else
    upperLimit = maximumValue + max(0.005, 0.05 * maximumValue);
end
limits = [0, upperLimit];
end

function plot_failed_magnitudes( ...
    ax, paired, methodName, sourceOffset, methodColour)
%PLOT_FAILED_MAGNITUDES Mark finite searches that failed direct replay.

switch methodName
    case 'Shinbrot'
        failed = ~logical(paired.ShinbrotReplayHit) & ...
            isfinite(paired.ShinbrotAbsoluteP);
        values = paired.ShinbrotAbsoluteP;
        label = 'Shinbrot replay failed';
    case 'Hybrid'
        failed = ~logical(paired.HybridReplayHit) & ...
            isfinite(paired.HybridAbsoluteP);
        values = paired.HybridAbsoluteP;
        label = 'Hybrid replay failed';
    otherwise
        error('Unknown method name: %s', methodName);
end

if any(failed)
    plot(ax, paired.SourceIndex(failed) + sourceOffset, values(failed), ...
        'x', 'LineStyle', 'none', 'Color', methodColour, ...
        'LineWidth', 1.5, 'MarkerSize', 7.5, ...
        'DisplayName', label);
end
end

function show_legend(ax, location, orientation)
%SHOW_LEGEND Display only plotted series with non-empty legend labels.

objects = findobj(ax, '-property', 'DisplayName');
labels = get(objects, 'DisplayName');
if ischar(labels)
    labels = {labels};
end
if isempty(labels) || ~any(~cellfun('isempty', labels))
    return;
end
legend(ax, 'show', 'Location', location, ...
    'Orientation', orientation, 'Box', 'off', ...
    'FontSize', 10);
end

function export_ensemble_svgs(figures, outputDirectory)
%EXPORT_ENSEMBLE_SVGS Save all figures as editable vector graphics.

if isstring(outputDirectory)
    if ~isscalar(outputDirectory)
        error('outputDirectory must be a character vector or scalar string.');
    end
    outputDirectory = char(outputDirectory);
end
if ~ischar(outputDirectory)
    error('outputDirectory must be a character vector or scalar string.');
end
if ~exist(outputDirectory, 'dir')
    [created, message] = mkdir(outputDirectory);
    if ~created
        error('Could not create output directory: %s', message);
    end
end

fileNames = { ...
    '01_returned_perturbation_by_source.svg', ...
    '02_paired_perturbation_difference.svg', ...
    '03_selected_targeting_horizon.svg', ...
    '04_paired_perturbation_comparison.svg'};
for figureNumber = 1:numel(figures)
    outputFile = fullfile(outputDirectory, fileNames{figureNumber});
    print(figures(figureNumber), outputFile, '-dsvg');
end
end
