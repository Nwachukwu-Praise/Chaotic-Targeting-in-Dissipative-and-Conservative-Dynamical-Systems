function figureHandles = so_plot_results(result)
%SO_PLOT_RESULTS Figures for the corrected Schroer-Ott demonstration.
cfg = result.configuration;
yWindow = visible_y_window(cfg);
figureHandles = gobjects(0);

fig1 = figure('Name', 'Schroer-Ott phase portrait and route', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 11]);
ax = axes(fig1);
hold(ax, 'on');
plot_phase_portrait_background(ax, cfg, yWindow);
plot_rectangle(ax, cfg.sourceRectangle, [0.1 0.35 0.9], 'source');
plot_rectangle(ax, cfg.targetRectangle, [0.9 0.15 0.15], 'target');
plot_orbit_catalogue(ax, result.orbitCatalogue.chains, [0.45 0.45 0.45]);
plot_route_chains(ax, result.route.chains);
plot_proxy_components(ax, result.route);
axis(ax, [0 1 yWindow]);
xlabel(ax, 'x mod 1');
ylabel(ax, 'y');
title(ax, 'Uncontrolled phase portrait, retained UPO chains, and operational route');
grid(ax, 'on');
legend(ax, 'Location', 'eastoutside');
hold(ax, 'off');
save_figure(fig1, cfg, 'phase_portrait_route');
figureHandles(end + 1) = fig1; %#ok<AGROW>

fig2 = figure('Name', 'Independently replayed controlled trajectory', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 11]);
layout = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axMain = nexttile(layout, 1);
hold(axMain, 'on');
plot_phase_portrait_background(axMain, cfg, yWindow);
plot_rectangle(axMain, cfg.sourceRectangle, [0.1 0.35 0.9], 'source');
plot_rectangle(axMain, cfg.targetRectangle, [0.9 0.15 0.15], 'target');
plot_proxy_components(axMain, result.route);
plot_executed_segments(axMain, result.executionSegments, cfg);
axis(axMain, [0 1 yWindow]);
xlabel(axMain, 'x mod 1');
ylabel(axMain, 'y');
title(axMain, 'Replay trajectory with instantaneous y-kicks');
grid(axMain, 'on');
legend(axMain, 'Location', 'best');

axInset = nexttile(layout, 2);
hold(axInset, 'on');
plot_rectangle(axInset, cfg.targetRectangle, [0.9 0.15 0.15], 'target');
plot_final_entry(axInset, result.executionSegments(end), cfg);
pad = 0.025;
axis(axInset, [cfg.targetRectangle.xMin - pad, cfg.targetRectangle.xMax + pad, ...
    cfg.targetRectangle.yMin - pad, cfg.targetRectangle.yMax + pad]);
xlabel(axInset, 'x mod 1');
ylabel(axInset, 'y');
title(axInset, sprintf('Final containment: %d', result.targetContained));
grid(axInset, 'on');
hold(axInset, 'off');
save_figure(fig2, cfg, 'controlled_trajectory_replay');
figureHandles(end + 1) = fig2; %#ok<AGROW>

fig3 = figure('Name', 'Switch objectives', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 8.5]);
stageCount = max(1, numel(result.switchEvaluations));
tiledlayout(fig3, stageCount, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(result.switchEvaluations)
    ax = nexttile;
    probes = result.switchEvaluations(i).probes;
    j = [probes.j];
    J = [probes.objectiveJ];
    finite = [probes.finite];
    plot(ax, j(finite), J(finite), 'ko-', 'MarkerSize', 4, ...
        'DisplayName', 'finite probe');
    hold(ax, 'on');
    if any(~finite)
        plot(ax, j(~finite), zeros(1, sum(~finite)), 'rx', 'DisplayName', 'unresolved/no route');
    end
    selected = result.switchEvaluations(i).selected;
    if selected.finite
        plot(ax, selected.j, selected.objectiveJ, 'o', 'Color', [0.9 0.15 0.15], ...
            'MarkerFaceColor', [0.9 0.15 0.15], 'DisplayName', 'selected');
    end
    ylabel(ax, sprintf('stage %d', result.switchEvaluations(i).stage));
    grid(ax, 'on');
    legend(ax, 'Location', 'best');
end
xlabel('candidate switch index j');
save_figure(fig3, cfg, 'switch_objectives');
figureHandles(end + 1) = fig3; %#ok<AGROW>

if isfield(result, 'diagnosticManifolds') && ~isempty(result.diagnosticManifolds)
    fig4 = figure('Name', 'Diagnostic manifolds and replayed switches', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 17 11]);
    ax = axes(fig4);
    hold(ax, 'on');
    plot_phase_portrait_background(ax, cfg, yWindow);
    plot_diagnostic_manifolds(ax, result.diagnosticManifolds);
    plot_rectangle(ax, cfg.sourceRectangle, [0.1 0.35 0.9], 'source');
    plot_rectangle(ax, cfg.targetRectangle, [0.9 0.15 0.15], 'target');
    plot_executed_segments(ax, result.executionSegments, cfg);
    axis(ax, [0 1 yWindow]);
    xlabel(ax, 'x mod 1');
    ylabel(ax, 'y');
    title(ax, 'Stable/unstable manifold diagnostics with replayed switches');
    grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside');
    hold(ax, 'off');
    save_figure(fig4, cfg, 'diagnostic_manifolds_switch_points');
    figureHandles(end + 1) = fig4; %#ok<AGROW>
end
end

function yWindow = visible_y_window(cfg)
yMin = min(cfg.sourceRectangle.yMin, cfg.targetRectangle.yMin) - 0.12;
yMax = max(cfg.sourceRectangle.yMax, cfg.targetRectangle.yMax) + 0.12;
yWindow = [yMin, yMax];
end

function plot_phase_portrait_background(ax, cfg, yWindow)
xSeeds = linspace(0, 1, 37);
xSeeds(end) = [];
ySeeds = linspace(yWindow(1), yWindow(2), 13);
points = zeros(2, numel(xSeeds) * numel(ySeeds) * 80);
idx = 0;
for ix = 1:numel(xSeeds)
    for iy = 1:numel(ySeeds)
        z = [xSeeds(ix); ySeeds(iy)];
        for n = 1:80
            z = so_standard_map_lifted(z, cfg);
            if z(2) >= yWindow(1) && z(2) <= yWindow(2)
                idx = idx + 1;
                points(:, idx) = [so_wrap_x(z(1)); z(2)];
            end
        end
    end
end
points = points(:, 1:idx);
scatter(ax, points(1, :), points(2, :), 2, [0.78 0.78 0.78], '.', ...
    'DisplayName', 'uncontrolled iterates');
end

function plot_orbit_catalogue(ax, chains, color)
for i = 1:numel(chains)
    pts = chains(i).pointsCylinder;
    if strcmp(chains(i).classification, 'direct-hyperbolic')
        marker = 'x';
    elseif strcmp(chains(i).classification, 'inverse-hyperbolic')
        marker = '+';
    else
        marker = '.';
    end
    plot(ax, pts(1, :), pts(2, :), marker, 'Color', color, ...
        'MarkerSize', 5, 'HandleVisibility', 'off');
end
end

function plot_route_chains(ax, chains)
for i = 1:numel(chains)
    pts = chains(i).pointsCylinder;
    plot(ax, pts(1, :), pts(2, :), 'ko', 'MarkerSize', 5, ...
        'LineWidth', 1.1, 'DisplayName', sprintf('\\omega=%.3g', chains(i).omega));
    text(ax, pts(1, 1), pts(2, 1), sprintf(' %.3g', chains(i).omega), 'FontSize', 8);
end
end

function plot_proxy_components(ax, route)
for i = 1:numel(route.targetComponents)
    comps = route.targetComponents{i};
    for c = 1:numel(comps)
        s = linspace(0, 1, 80);
        pts = comps{c}.gamma0(s);
        pts(1, :) = so_wrap_x(pts(1, :));
        plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.35 0.35 0.35], ...
            'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
end
end

function plot_executed_segments(ax, segments, cfg)
for i = 1:numel(segments)
    seg = segments(i);
    pre = so_to_cylinder(seg.preControlState);
    post = so_to_cylinder(seg.postControlState);
    quiver(ax, pre(1), pre(2), 0, post(2) - pre(2), 0, ...
        'Color', [0.0 0.45 0.75], 'LineWidth', 1.1, 'MaxHeadSize', 1.4, ...
        'DisplayName', control_label(i, 'control kick'));
    plot(ax, pre(1), pre(2), 's', 'Color', [0.0 0.45 0.75], ...
        'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'DisplayName', control_label(i, 'pre-control'));
    plot(ax, post(1), post(2), 'o', 'Color', [0.0 0.45 0.75], ...
        'MarkerFaceColor', [0.0 0.45 0.75], 'MarkerSize', 4, ...
        'DisplayName', control_label(i, 'post-control'));
    pts = break_periodic_line(seg.path);
    plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.05 0.05 0.05], ...
        'LineWidth', 1.2, 'DisplayName', control_label(i, 'replayed iterates'));
    sw = so_to_cylinder(seg.switchState);
    plot(ax, sw(1), sw(2), 'd', 'Color', [0.85 0.2 0.0], ...
        'MarkerFaceColor', [0.85 0.2 0.0], 'MarkerSize', 5, ...
        'DisplayName', control_label(i, 'switch/final state'));
end
finalState = so_to_cylinder(segments(end).finalState);
plot(ax, finalState(1), finalState(2), 'p', 'Color', [0.5 0 0.6], ...
    'MarkerFaceColor', [0.5 0 0.6], 'MarkerSize', 9, 'DisplayName', 'final state');
end

function label = control_label(i, textLabel)
if i == 1
    label = textLabel;
else
    label = '';
end
end

function plot_final_entry(ax, segment, cfg)
pts = break_periodic_line(segment.path);
plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.05 0.05 0.05], ...
    'LineWidth', 1.3, 'DisplayName', 'final replay segment');
finalState = so_to_cylinder(segment.finalState);
plot(ax, finalState(1), finalState(2), 'p', 'Color', [0.5 0 0.6], ...
    'MarkerFaceColor', [0.5 0 0.6], 'MarkerSize', 9, 'DisplayName', 'final state');
legend(ax, 'Location', 'best');
end

function pts = break_periodic_line(path)
pts = so_to_cylinder(path);
if size(pts, 2) < 2
    return;
end
out = pts(:, 1);
for i = 2:size(pts, 2)
    if abs(pts(1, i) - pts(1, i - 1)) > 0.5
        out = [out, [NaN; NaN]]; %#ok<AGROW>
    end
    out = [out, pts(:, i)]; %#ok<AGROW>
end
pts = out;
end

function plot_diagnostic_manifolds(ax, manifolds)
for i = 1:numel(manifolds)
    curve = manifolds(i).curve;
    pts = break_periodic_line(curve.pointsLifted);
    if manifolds(i).branchType == "unstable"
        col = [0.8 0.1 0.1];
        label = 'unstable manifold sample';
    else
        col = [0.1 0.2 0.8];
        label = 'stable manifold sample';
    end
    if i > 2
        label = '';
    end
    plot(ax, pts(1, :), pts(2, :), '-', 'Color', col, 'LineWidth', 0.7, ...
        'DisplayName', label);
end
end

function plot_rectangle(ax, rect, color, label)
if rect.xMax >= rect.xMin
    xs = [rect.xMin rect.xMax rect.xMax rect.xMin rect.xMin];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
else
    xs = [rect.xMin 1 1 rect.xMin rect.xMin NaN 0 rect.xMax rect.xMax 0 0];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin NaN ...
        rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
end
plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 1.8, 'DisplayName', label);
end

function save_figure(fig, cfg, baseName)
if ~cfg.saveFigures
    return;
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end
pngPath = fullfile(cfg.figureDirectory, [char(baseName), '.png']);
svgPath = fullfile(cfg.figureDirectory, [char(baseName), '.svg']);
epsPath = fullfile(cfg.figureDirectory, [char(baseName), '.eps']);
exportgraphics(fig, pngPath, 'Resolution', 180);
exportgraphics(fig, svgPath, 'ContentType', 'vector');
print(fig, epsPath, '-depsc', '-painters');
end
