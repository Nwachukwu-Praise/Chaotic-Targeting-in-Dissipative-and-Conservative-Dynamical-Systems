function fig = so_plot_ensemble(ensemble, saveFigure)
%SO_PLOT_ENSEMBLE Spread of controlled transfer times over random sources.
%
% Four panels:
%   (a) where the source points were drawn, on the phase portrait, coloured
%       by controlled iteration count -- so a systematic dependence on
%       position inside the source region would be visible immediately;
%   (b) histogram of controlled iteration counts;
%   (c) controlled iterations against source y, since the rotation bracket
%       (and therefore the number of resonances crossed) depends on it;
%   (d) controlled against uncontrolled transport time on a log axis.
if nargin < 2
    saveFigure = true;
end
cfg = ensemble.configuration;
trials = ensemble.trials;
ok = trials.success & isfinite(trials.totalControlledIterations);

fig = figure('Name', 'Source ensemble', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1 1 24 17]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) source positions ------------------------------------------------
ax1 = nexttile(layout, 1);
hold(ax1, 'on');
yWindow = so_plot_y_window(cfg);
so_plot_phase_background(ax1, cfg, yWindow, false);
draw_rect(ax1, cfg.sourceRectangle, [0.00 0.35 0.85]);
draw_rect(ax1, cfg.targetRectangle, [0.85 0.10 0.10]);
if any(ok)
    scatter(ax1, trials.sourceX(ok), trials.sourceY(ok), 34, ...
        trials.totalControlledIterations(ok), 'filled', 'MarkerEdgeColor', 'k');
    cb = colorbar(ax1);
    cb.Label.String = 'controlled iterations';
end
if any(~ok)
    plot(ax1, trials.sourceX(~ok), trials.sourceY(~ok), 'rx', ...
        'MarkerSize', 9, 'LineWidth', 1.4);
end
axis(ax1, [0 1 yWindow]);
xlabel(ax1, 'x  (mod 1)'); ylabel(ax1, 'y');
title(ax1, sprintf('(a) %d source points  (%d reached target)', ...
    height(trials), sum(ok)));
grid(ax1, 'on'); hold(ax1, 'off');

% ---- (b) histogram -------------------------------------------------------
ax2 = nexttile(layout, 2);
if any(ok)
    counts = trials.totalControlledIterations(ok);
    edges = (min(counts) - 0.5):1:(max(counts) + 0.5);
    if numel(edges) < 2
        edges = [min(counts) - 0.5, min(counts) + 0.5];
    end
    histogram(ax2, counts, edges, 'FaceColor', [0.00 0.35 0.85], ...
        'EdgeColor', 'w');
    hold(ax2, 'on');
    xline(ax2, median(counts), '-', sprintf('median %g', median(counts)), ...
        'Color', [0.85 0.35 0.00], 'LineWidth', 1.6, ...
        'LabelVerticalAlignment', 'top');
    hold(ax2, 'off');
end
xlabel(ax2, 'controlled iterations to target');
ylabel(ax2, 'number of source points');
title(ax2, '(b) spread of controlled transfer time');
grid(ax2, 'on');

% ---- (c) dependence on source y -----------------------------------------
ax3 = nexttile(layout, 3);
hold(ax3, 'on');
if any(ok)
    gscatter_by_route(ax3, trials(ok, :));
end
if any(~ok)
    yl = ylim(ax3);
    plot(ax3, trials.sourceY(~ok), repmat(yl(1), sum(~ok), 1), 'rx', ...
        'MarkerSize', 9, 'LineWidth', 1.4, 'DisplayName', 'not reached');
end
xlabel(ax3, 'source y');
ylabel(ax3, 'controlled iterations');
title(ax3, '(c) transfer time against source y, by route');
grid(ax3, 'on');
legend(ax3, 'Location', 'best');
hold(ax3, 'off');

% ---- (d) controlled versus uncontrolled ---------------------------------
ax4 = nexttile(layout, 4);
haveUnc = ok & isfinite(trials.uncontrolledTransportTime);
if any(haveUnc)
    semilogy(ax4, trials.totalControlledIterations(haveUnc), ...
        trials.uncontrolledTransportTime(haveUnc), 'o', ...
        'Color', [0.00 0.35 0.85], 'MarkerFaceColor', [0.60 0.75 0.95], ...
        'DisplayName', 'source point');
    hold(ax4, 'on');
    xl = xlim(ax4);
    plot(ax4, xl, xl, 'k--', 'DisplayName', 'no gain');
    hold(ax4, 'off');
    legend(ax4, 'Location', 'best');
    ylabel(ax4, 'uncontrolled first passage (iterations)');
else
    text(0.5, 0.5, 'uncontrolled times not measured', ...
        'Parent', ax4, 'HorizontalAlignment', 'center');
end
xlabel(ax4, 'controlled iterations');
title(ax4, '(d) what the control buys');
grid(ax4, 'on');

title(layout, sprintf('Source ensemble, %s case, k = %g, |\\delta_n| \\leq %g', ...
    char(cfg.caseName), cfg.k, cfg.controlAmplitude));

if saveFigure
    exportCfg = cfg;
    exportCfg.figureDirectory = fullfile(cfg.outputDirectory, 'ensemble', 'figures');
    so_export_figure(fig, exportCfg, 'source_ensemble');
end
end

function gscatter_by_route(ax, t)
routes = unique(t.routeOmegas);
colors = lines(max(1, numel(routes)));
for i = 1:numel(routes)
    m = t.routeOmegas == routes(i);
    plot(ax, t.sourceY(m), t.totalControlledIterations(m), 'o', ...
        'Color', colors(i, :), 'MarkerFaceColor', colors(i, :), ...
        'MarkerSize', 6, 'DisplayName', char(routes(i)));
end
end

function draw_rect(ax, rect, color)
if rect.xMax >= rect.xMin
    xs = [rect.xMin rect.xMax rect.xMax rect.xMin rect.xMin];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
else
    xs = [rect.xMin 1 1 rect.xMin rect.xMin NaN 0 rect.xMax rect.xMax 0 0];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin NaN ...
        rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
end
plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 1.8, 'HandleVisibility', 'off');
end
