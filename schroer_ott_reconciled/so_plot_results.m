function figureHandles = so_plot_results(result)
%SO_PLOT_RESULTS Figures for one targeting case, over the phase portrait.
%
% Every phase-space panel is drawn on top of the uncontrolled standard map
% (so_plot_phase_background), so the route, the proxy balls and the
% controlled trajectory can be read against the island/chaos structure they
% exploit rather than against an empty axis.
cfg = result.configuration;
yWindow = so_plot_y_window(cfg);
figureHandles = gobjects(0);

% ---------------------------------------------------------------- figure 1
fig1 = figure('Name', 'Phase portrait and route', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 19 12]);
ax = axes(fig1); %#ok<LAXES>
hold(ax, 'on');
so_plot_phase_background(ax, cfg, yWindow, true);
plot_orbit_catalogue(ax, result.orbitCatalogue.chains);
plot_route_chains(ax, result.route);
plot_proxy_components(ax, result.route, cfg);
plot_rectangle(ax, cfg.sourceRectangle, [0.00 0.35 0.85], 'source region');
plot_rectangle(ax, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
axis(ax, [0 1 yWindow]);
xlabel(ax, 'x  (mod 1)');
ylabel(ax, 'y');
title(ax, sprintf('%s  --  route %s', char(cfg.caseLabel), ...
    route_label(result.route)));
grid(ax, 'on');
legend(ax, 'Location', 'eastoutside');
hold(ax, 'off');
save_figure(fig1, cfg, 'phase_portrait_route');
figureHandles(end + 1) = fig1;

% ---------------------------------------------------------------- figure 2
fig2 = figure('Name', 'Controlled trajectory', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 21 11]);
layout = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
axMain = nexttile(layout, 1);
hold(axMain, 'on');
so_plot_phase_background(axMain, cfg, yWindow, true);
plot_proxy_components(axMain, result.route, cfg);
plot_rectangle(axMain, cfg.sourceRectangle, [0.00 0.35 0.85], 'source region');
plot_rectangle(axMain, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
plot_execution(axMain, result);
axis(axMain, [0 1 yWindow]);
xlabel(axMain, 'x  (mod 1)');
ylabel(axMain, 'y');
title(axMain, sprintf('%d controls, %d map iterations', ...
    result.numberOfControls, result.totalExecutedIterations));
grid(axMain, 'on');
legend(axMain, 'Location', 'southoutside', 'NumColumns', 2);

axInset = nexttile(layout, 2);
hold(axInset, 'on');
plot_target_closeup(axInset, result, cfg);
title(axInset, sprintf('Final approach  (contained: %d)', result.targetContained));
grid(axInset, 'on');
hold(axInset, 'off');
save_figure(fig2, cfg, 'controlled_trajectory');
figureHandles(end + 1) = fig2;

% ---------------------------------------------------------------- figure 3
if ~isempty(result.switchEvaluations)
    fig3 = figure('Name', 'Switch objectives', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 17 4 + 3 * numel(result.switchEvaluations)]);
    tiledlayout(fig3, numel(result.switchEvaluations), 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(result.switchEvaluations)
        axJ = nexttile;
        plot_switch_objective(axJ, result.switchEvaluations(i));
    end
    xlabel('candidate switch index  j');
    save_figure(fig3, cfg, 'switch_objectives');
    figureHandles(end + 1) = fig3;
end

% ---------------------------------------------------------------- figure 4
if isfield(result, 'diagnosticManifolds') && ~isempty(result.diagnosticManifolds)
    fig4 = figure('Name', 'Manifold diagnostics', 'Color', 'w', ...
        'Units', 'centimeters', 'Position', [2 2 19 12]);
    ax4 = axes(fig4); %#ok<LAXES>
    hold(ax4, 'on');
    so_plot_phase_background(ax4, cfg, yWindow, true);
    plot_diagnostic_manifolds(ax4, result.diagnosticManifolds);
    plot_rectangle(ax4, cfg.sourceRectangle, [0.00 0.35 0.85], 'source region');
    plot_rectangle(ax4, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
    plot_execution(ax4, result);
    axis(ax4, [0 1 yWindow]);
    xlabel(ax4, 'x  (mod 1)');
    ylabel(ax4, 'y');
    title(ax4, 'A posteriori manifold branches and executed switch points');
    grid(ax4, 'on');
    legend(ax4, 'Location', 'eastoutside');
    hold(ax4, 'off');
    save_figure(fig4, cfg, 'diagnostic_manifolds_switch_points');
    figureHandles(end + 1) = fig4;
end
end

% ======================================================================= %

function label = route_label(route)
if isempty(route.rotationNumbers)
    label = '(empty rotation bracket: single forward-backward step)';
else
    parts = arrayfun(@(w) sprintf('%.4g', w), route.rotationNumbers, ...
        'UniformOutput', false);
    label = sprintf('omega = %s', strjoin(parts, ', '));
end
end

function plot_orbit_catalogue(ax, chains)
shown = false(1, 3);
for i = 1:numel(chains)
    pts = chains(i).pointsCylinder;
    switch chains(i).classification
        case 'direct-hyperbolic'
            marker = 'x'; col = [0.10 0.10 0.10]; slot = 1;
            name = 'direct-hyperbolic UPO';
        case 'inverse-hyperbolic'
            marker = '+'; col = [0.45 0.20 0.55]; slot = 2;
            name = 'inverse-hyperbolic UPO';
        otherwise
            marker = 'o'; col = [0.20 0.45 0.20]; slot = 3;
            name = 'elliptic centre';
    end
    if shown(slot)
        args = {'HandleVisibility', 'off'};
    else
        args = {'DisplayName', name};
        shown(slot) = true;
    end
    plot(ax, pts(1, :), pts(2, :), marker, 'Color', col, ...
        'MarkerSize', 5, 'LineWidth', 0.9, args{:});
end
end

function plot_route_chains(ax, route)
for i = 1:numel(route.chains)
    pts = route.chains(i).pointsCylinder;
    plot(ax, pts(1, :), pts(2, :), 'x', 'Color', [0.85 0.35 0.00], ...
        'MarkerSize', 10, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('route \\omega = %.4g', route.chains(i).omega));
end
end

function plot_proxy_components(ax, route, cfg) %#ok<INUSD>
first = true;
for i = 1:numel(route.targetComponents)
    comps = route.targetComponents{i};
    for c = 1:numel(comps)
        s = linspace(0, 1, 121);
        pts = comps{c}.gamma0(s);
        pts(1, :) = so_wrap_x(pts(1, :));
        pts = break_periodic_line(pts);
        if first
            args = {'DisplayName', 'proxy target ball'};
            first = false;
        else
            args = {'HandleVisibility', 'off'};
        end
        plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.85 0.35 0.00], ...
            'LineWidth', 0.8, args{:});
    end
end
end

function plot_execution(ax, result)
%PLOT_EXECUTION Iterates as solid dots, controls as open circles.
%
% This follows the paper's own caption for Figure 3: "The solid circles
% denote the trajectory of length 125 with controls applied at open
% circles."  Successive iterates of the standard map are far apart -- x
% advances by roughly y each step -- so joining them with line segments
% draws chords across the whole cell and buries the structure underneath in
% spaghetti.  The iterates themselves are the trajectory.
segments = result.executionSegments;
if isempty(segments)
    z = so_to_cylinder(result.sourceState);
    plot(ax, z(1), z(2), 's', 'Color', [0.00 0.35 0.85], 'MarkerFaceColor', 'w', ...
        'MarkerSize', 8, 'LineWidth', 1.3, 'DisplayName', 'source state');
    return;
end

pts = zeros(2, 0);
ctrl = zeros(2, 0);
for i = 1:numel(segments)
    pts = [pts, so_to_cylinder(segments(i).path)]; %#ok<AGROW>
    ctrl = [ctrl, so_to_cylinder(segments(i).postControlState)]; %#ok<AGROW>
end
plot(ax, pts(1, :), pts(2, :), 'o', 'Color', [0.05 0.05 0.05], ...
    'MarkerFaceColor', [0.05 0.05 0.05], 'MarkerSize', 3.5, 'LineStyle', 'none', ...
    'DisplayName', 'controlled iterates');
plot(ax, ctrl(1, :), ctrl(2, :), 'o', 'Color', [0.00 0.45 0.75], ...
    'MarkerFaceColor', 'none', 'MarkerSize', 9, 'LineWidth', 1.5, ...
    'LineStyle', 'none', 'DisplayName', 'control applied');
z0 = so_to_cylinder(result.sourceState);
plot(ax, z0(1), z0(2), 's', 'Color', [0.00 0.35 0.85], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, 'LineWidth', 1.3, 'DisplayName', 'source state');
zf = so_to_cylinder(result.finalState);
plot(ax, zf(1), zf(2), 'p', 'Color', [0.45 0.00 0.55], ...
    'MarkerFaceColor', [0.45 0.00 0.55], 'MarkerSize', 11, ...
    'DisplayName', 'final state');
end

function plot_target_closeup(ax, result, cfg)
plot_rectangle(ax, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
segments = result.executionSegments;
if ~isempty(segments)
    seg = segments(end);
    pts = so_to_cylinder(seg.path);
    plot(ax, pts(1, :), pts(2, :), 'o', 'Color', [0.05 0.05 0.05], ...
        'MarkerFaceColor', [0.05 0.05 0.05], 'MarkerSize', 4.5, ...
        'LineStyle', 'none', 'DisplayName', 'final segment iterates');
    post = so_to_cylinder(seg.postControlState);
    plot(ax, post(1), post(2), 'o', 'Color', [0.00 0.45 0.75], ...
        'MarkerFaceColor', 'none', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'LineStyle', 'none', 'DisplayName', 'control applied');
end
zf = so_to_cylinder(result.finalState);
plot(ax, zf(1), zf(2), 'p', 'Color', [0.45 0.00 0.55], ...
    'MarkerFaceColor', [0.45 0.00 0.55], 'MarkerSize', 12, ...
    'DisplayName', 'final state');
pad = 0.5 * max(cfg.targetRectangle.xMax - cfg.targetRectangle.xMin, ...
    cfg.targetRectangle.yMax - cfg.targetRectangle.yMin);
axis(ax, [cfg.targetRectangle.xMin - pad, cfg.targetRectangle.xMax + pad, ...
    cfg.targetRectangle.yMin - pad, cfg.targetRectangle.yMax + pad]);
xlabel(ax, 'x  (mod 1)');
ylabel(ax, 'y');
legend(ax, 'Location', 'best');
end

function plot_switch_objective(ax, evaluation)
probes = evaluation.probes;
j = [probes.j];
J = [probes.objectiveJ];
finite = [probes.finite];
pruned = [probes.pruned];
hold(ax, 'on');
if any(finite)
    plot(ax, j(finite), J(finite), 'o-', 'Color', [0.10 0.10 0.10], ...
        'MarkerSize', 4, 'MarkerFaceColor', 'w', 'DisplayName', 'J(j) = j + \tau_{res}');
end
unres = ~finite & ~pruned;
if any(unres)
    yl = ylim(ax);
    plot(ax, j(unres), repmat(yl(1), 1, sum(unres)), 'x', ...
        'Color', [0.85 0.10 0.10], 'DisplayName', 'no resolved connection');
end
if any(pruned) && any(finite)
    plot(ax, j(pruned), repmat(min(J(finite)), 1, sum(pruned)), '.', ...
        'Color', [0.60 0.60 0.60], 'DisplayName', 'pruned (j >= J_{best})');
end
sel = evaluation.selected;
if sel.finite
    plot(ax, sel.j, sel.objectiveJ, 'o', 'Color', [0.85 0.35 0.00], ...
        'MarkerFaceColor', [0.85 0.35 0.00], 'MarkerSize', 7, ...
        'DisplayName', 'selected switch');
end
ylabel(ax, sprintf('stage %d', evaluation.stage));
title(ax, sprintf('%s  \\rightarrow  %s', char(evaluation.currentTargetID), ...
    char(evaluation.nextTargetID)), 'FontWeight', 'normal');
grid(ax, 'on');
legend(ax, 'Location', 'best');
hold(ax, 'off');
end

function plot_diagnostic_manifolds(ax, manifolds)
shownU = false;
shownS = false;
for i = 1:numel(manifolds)
    pts = break_periodic_line(manifolds(i).curve.pointsLifted);
    if manifolds(i).branchType == "unstable"
        col = [0.80 0.10 0.10];
        if shownU
            args = {'HandleVisibility', 'off'};
        else
            args = {'DisplayName', 'unstable manifold'};
            shownU = true;
        end
    else
        col = [0.10 0.20 0.80];
        if shownS
            args = {'HandleVisibility', 'off'};
        else
            args = {'DisplayName', 'stable manifold'};
            shownS = true;
        end
    end
    plot(ax, pts(1, :), pts(2, :), '-', 'Color', col, 'LineWidth', 0.7, args{:});
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
plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 2.0, 'DisplayName', label);
end

function pts = break_periodic_line(path)
%BREAK_PERIODIC_LINE Insert NaN where a line wraps the x-period.
pts = so_to_cylinder(path);
if size(pts, 2) < 2
    return;
end
jump = [false, abs(diff(pts(1, :))) > 0.5];
out = zeros(2, size(pts, 2) + sum(jump));
k = 0;
for i = 1:size(pts, 2)
    if jump(i)
        k = k + 1;
        out(:, k) = [NaN; NaN];
    end
    k = k + 1;
    out(:, k) = pts(:, i);
end
pts = out(:, 1:k);
end

function save_figure(fig, cfg, baseName)
if ~cfg.saveFigures
    return;
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end
so_export_figure(fig, cfg, baseName);
end
