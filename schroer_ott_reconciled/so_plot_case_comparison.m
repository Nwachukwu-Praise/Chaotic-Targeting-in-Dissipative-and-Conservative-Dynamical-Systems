function fig = so_plot_case_comparison(results, saveFigure)
%SO_PLOT_CASE_COMPARISON The three case studies on one sheet.
%
% Each panel carries the same uncontrolled phase portrait underneath, so the
% three transfers can be compared against identical structure: the reader
% can see that the horizontal case never has to cross a resonance, while the
% vertical case crosses three.
if nargin < 2
    saveFigure = true;
end
n = numel(results);
fig = figure('Name', 'Case comparison', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1 1 7 * n + 3, 11]);
layout = tiledlayout(fig, 1, n, 'TileSpacing', 'compact', 'Padding', 'compact');

axList = gobjects(1, n);
for i = 1:n
    r = results{i};
    cfg = r.configuration;
    ax = nexttile(layout, i);
    axList(i) = ax;
    hold(ax, 'on');
    yWindow = so_plot_y_window(cfg);
    so_plot_phase_background(ax, cfg, yWindow, i == 1);
    draw_rect(ax, cfg.sourceRectangle, [0.00 0.35 0.85], 'source region', i == 1);
    draw_rect(ax, cfg.targetRectangle, [0.85 0.10 0.10], 'target region', i == 1);
    draw_proxies(ax, r.route, i == 1);
    draw_path(ax, r, i == 1);
    axis(ax, [0 1 yWindow]);
    xlabel(ax, 'x  (mod 1)');
    if i == 1
        ylabel(ax, 'y');
    end
    if isempty(r.route.rotationNumbers)
        routeText = 'no resonance crossed';
    else
        routeText = sprintf('%d resonance(s) crossed', numel(r.route.rotationNumbers));
    end
    title(ax, {upper(char(cfg.caseName)), ...
        sprintf('%s, %d controls, %d iterations', routeText, ...
        r.numberOfControls, r.totalExecutedIterations)}, 'FontWeight', 'normal');
    grid(ax, 'on');
    hold(ax, 'off');
end
lgd = legend(axList(1), 'Orientation', 'horizontal', 'NumColumns', 3);
lgd.Layout.Tile = 'south';
title(layout, 'Schroer-Ott pass targeting: three source/target geometries, k = 1.25');

if saveFigure
    cfg = results{1}.configuration;
    exportCfg = cfg;
    exportCfg.figureDirectory = fullfile(pwd, 'outputs', 'figures');
    so_export_figure(fig, exportCfg, 'case_comparison');
end
end

function draw_rect(ax, rect, color, label, showLegend)
if rect.xMax >= rect.xMin
    xs = [rect.xMin rect.xMax rect.xMax rect.xMin rect.xMin];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
else
    xs = [rect.xMin 1 1 rect.xMin rect.xMin NaN 0 rect.xMax rect.xMax 0 0];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin NaN ...
        rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
end
args = legend_args(label, showLegend);
plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 2.0, args{:});
end

function draw_proxies(ax, route, showLegend)
first = true;
for i = 1:numel(route.targetComponents)
    comps = route.targetComponents{i};
    for c = 1:numel(comps)
        pts = comps{c}.gamma0(linspace(0, 1, 121));
        pts(1, :) = so_wrap_x(pts(1, :));
        pts = split_wrap(pts);
        args = legend_args('proxy target ball', showLegend && first);
        plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.85 0.35 0.00], ...
            'LineWidth', 0.8, args{:});
        first = false;
    end
end
end

function draw_path(ax, r, showLegend)
first = true;
for i = 1:numel(r.executionSegments)
    seg = r.executionSegments(i);
    pts = split_wrap(so_to_cylinder(seg.path));
    args = legend_args('controlled trajectory', showLegend && first);
    plot(ax, pts(1, :), pts(2, :), '-', 'Color', [0.05 0.05 0.05], ...
        'LineWidth', 1.3, args{:});
    plot(ax, pts(1, :), pts(2, :), '.', 'Color', [0.05 0.05 0.05], ...
        'MarkerSize', 7, 'HandleVisibility', 'off');
    pre = so_to_cylinder(seg.preControlState);
    post = so_to_cylinder(seg.postControlState);
    args = legend_args('control kick', showLegend && first);
    quiver(ax, pre(1), pre(2), 0, post(2) - pre(2), 0, 'Color', [0.00 0.45 0.75], ...
        'LineWidth', 1.4, 'MaxHeadSize', 8, args{:});
    first = false;
end
z0 = so_to_cylinder(r.sourceState);
srcArgs = legend_args('source state', showLegend);
plot(ax, z0(1), z0(2), 's', 'Color', [0.00 0.35 0.85], 'MarkerFaceColor', 'w', ...
    'MarkerSize', 8, 'LineWidth', 1.3, srcArgs{:});
zf = so_to_cylinder(r.finalState);
finArgs = legend_args('final state', showLegend);
plot(ax, zf(1), zf(2), 'p', 'Color', [0.45 0.00 0.55], ...
    'MarkerFaceColor', [0.45 0.00 0.55], 'MarkerSize', 11, finArgs{:});
end

function args = legend_args(label, showLegend)
if showLegend
    args = {'DisplayName', label};
else
    args = {'HandleVisibility', 'off'};
end
end

function pts = split_wrap(pts)
if size(pts, 2) < 2
    return;
end
jump = [false, abs(diff(pts(1, :))) > 0.5];
out = zeros(2, size(pts, 2) + sum(jump));
k = 0;
for i = 1:size(pts, 2)
    if jump(i)
        k = k + 1; out(:, k) = [NaN; NaN];
    end
    k = k + 1; out(:, k) = pts(:, i);
end
pts = out(:, 1:k);
end
