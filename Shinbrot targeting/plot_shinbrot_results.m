function plot_shinbrot_results( ...
    mapData, uncontrolled, controlled, target, params, control)
%PLOT_SHINBROT_RESULTS Visualize the return map and targeting result.

figure('Name', 'Lorenz attractor and accepted section', 'Color', 'w');
plot3(mapData.X(:, 1), mapData.X(:, 2), mapData.X(:, 3), ...
    'Color', [0.08, 0.30, 0.68], 'LineWidth', 0.55);
hold on;
draw_section_half_plane(mapData.X, params);
plot3(target.state(1), target.state(2), target.state(3), ...
    'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 13);
grid on; box on; view(38, 24);
xlabel('x'); ylabel('y'); zlabel('z');
title(sprintf('Lorenz attractor and section z = %.3f, x > %.1f', ...
    params.zSection, params.xSectionMin));
legend({'Uncontrolled attractor', 'Accepted Poincare half-plane', ...
    'Target fixed point'}, 'Location', 'best');

figure('Name', 'Uncontrolled scalar return map', 'Color', 'w');
plot(mapData.xN, mapData.xNext, '.', ...
    'Color', [0.18, 0.43, 0.72], 'MarkerSize', 6);
hold on;
mapLimits = [min([mapData.xN, mapData.xNext]), ...
    max([mapData.xN, mapData.xNext])];
plot(mapLimits, mapLimits, 'k--', 'LineWidth', 1);
plot(target.x, target.x, 'rp', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 13);
grid on; box on; axis equal;
xlim(mapLimits); ylim(mapLimits);
xlabel('X_n');
ylabel('X_{n+1}');
title('Nominal Poincare return map, pControl = 0');
legend({'Return map', 'X_{n+1} = X_n', 'Target fixed point'}, ...
    'Location', 'best');

figure('Name', 'Controlled and uncontrolled section sequences', ...
    'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(0:(numel(uncontrolled.xSequence) - 1), ...
    uncontrolled.xSequence, '-o', ...
    'Color', [0.08, 0.35, 0.75], 'MarkerSize', 3, ...
    'MarkerFaceColor', [0.08, 0.35, 0.75]);
hold on;
draw_target_band(target, max(numel(uncontrolled.xSequence) - 1, 1));
xlabel('Accepted section-crossing number');
ylabel('X_n');
title(sprintf('Uncontrolled recurrence: %s crossings', ...
    format_count(uncontrolled.numCrossings)));
grid on; box on;

nexttile;
plot(0:(numel(controlled.xSequence) - 1), ...
    controlled.xSequence, '-o', ...
    'Color', [0.02, 0.55, 0.30], 'MarkerSize', 4, ...
    'MarkerFaceColor', [0.02, 0.55, 0.30]);
hold on;
draw_target_band(target, max(numel(controlled.xSequence) - 1, 1));
xlabel('Accepted section-crossing number');
ylabel('X_n');
title(sprintf('Controlled targeting: %s crossings', ...
    format_count(controlled.numCrossings)));
grid on; box on;

figure('Name', 'Controlled Lorenz trajectory', 'Color', 'w');
plot3(mapData.X(:, 1), mapData.X(:, 2), mapData.X(:, 3), ...
    'Color', [0.82, 0.82, 0.82], 'LineWidth', 0.4);
hold on;
if ~isempty(controlled.X)
    plot3(controlled.X(:, 1), controlled.X(:, 2), controlled.X(:, 3), ...
        'Color', [0.02, 0.55, 0.30], 'LineWidth', 1.4);
end
plot3(target.state(1), target.state(2), target.state(3), ...
    'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 13);
grid on; box on; view(38, 24);
xlabel('x'); ylabel('y'); zlabel('z');
title('Parameter-controlled approach to the figure-eight orbit region');
legend({'Lorenz attractor', 'Controlled trajectory', ...
    'Target fixed point'}, 'Location', 'best');

figure('Name', 'Parameter search', 'Color', 'w');
search = controlled.search;
hasGridCurve = any(isfinite(search.gridDistance));
hold on;
if hasGridCurve
    semilogy(search.pGrid, max(search.gridDistance, eps), '-o', ...
        'Color', [0.48, 0.16, 0.58], 'MarkerSize', 3, ...
        'MarkerFaceColor', [0.48, 0.16, 0.58]);
else
    text(0, 10 * target.tolerance, ...
        'A grid curve is unavailable for this full-interval method.', ...
        'HorizontalAlignment', 'center');
end
yline(target.tolerance, 'r--', 'Target tolerance', ...
    'LineWidth', 1.2);
if isfinite(search.bracketLow)
    xline(search.bracketLow, 'Color', [0.55, 0.55, 0.55], ...
        'LineStyle', ':', 'LineWidth', 1.1, 'Label', 'Bracket low');
end
if isfinite(search.bracketHigh)
    xline(search.bracketHigh, 'Color', [0.55, 0.55, 0.55], ...
        'LineStyle', ':', 'LineWidth', 1.1, 'Label', 'Bracket high');
end
if search.found
    xline(search.selectedP, 'k--', 'Selected pControl', ...
        'LineWidth', 1.2);
end
xlim([-control.deltaP, control.deltaP]);
lowerLimit = max(eps, target.tolerance / 20);
if hasGridCurve
    finiteDistance = search.gridDistance(isfinite(search.gridDistance));
    upperLimit = max([target.tolerance * 10, finiteDistance(:).']);
else
    upperLimit = target.tolerance * 100;
end
ylim([lowerLimit, upperLimit]);
xlabel('pControl');
ylabel('|X_n - X_t|');
title(sprintf('%s: %s, horizon n = %d', ...
    search.method, search.selectionMethod, search.plotHorizon));
grid on; box on;
end

function draw_section_half_plane(X, params)
%DRAW_SECTION_HALF_PLANE Show only the accepted x > xSectionMin portion.

xLimits = [params.xSectionMin, max(X(:, 1))];
yLimits = [min(X(:, 2)), max(X(:, 2))];
[xPlane, yPlane] = meshgrid(xLimits, yLimits);
zPlane = params.zSection * ones(size(xPlane));

surf(xPlane, yPlane, zPlane, ...
    'FaceColor', [0.95, 0.68, 0.10], ...
    'FaceAlpha', 0.22, 'EdgeColor', 'none');
end

function draw_target_band(target, maxCrossing)
%DRAW_TARGET_BAND Mark the one-dimensional target interval.

xPatch = [0, maxCrossing, maxCrossing, 0];
yPatch = target.x + target.tolerance * [-1, -1, 1, 1];
patch(xPatch, yPatch, [0.88, 0.12, 0.12], ...
    'FaceAlpha', 0.10, 'EdgeColor', [0.88, 0.12, 0.12], ...
    'LineWidth', 1);
yline(target.x, 'r:', 'X_t', 'LineWidth', 1);
end

function textValue = format_count(value)
%FORMAT_COUNT Format successful and failed hitting counts.

if isfinite(value)
    textValue = sprintf('%d', value);
else
    textValue = 'not hit';
end
end
