function fig = so_plot_noise(result, noise, saveFigure)
%SO_PLOT_NOISE Noisy replays against the deterministic trajectory.
%
%   (a) the headline noise level drawn over the phase portrait, every
%       realisation on top of the deterministic path, so the divergence can
%       be seen where it happens rather than only at the endpoint;
%   (b) close-up of the target region with every noisy endpoint;
%   (c) containment fraction against noise amplitude;
%   (d) endpoint displacement against noise amplitude, log-log.
if nargin < 3
    saveFigure = true;
end
cfg = result.configuration;
levels = noise.levels;

[~, headIdx] = min(abs(levels.sigmaOverDelta - noise.headlineFactor));
headRecord = noise.records(headIdx);

fig = figure('Name', 'Noise validation', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [1 1 24 17]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- (a) trajectories ----------------------------------------------------
ax1 = nexttile(layout, 1);
hold(ax1, 'on');
yWindow = so_plot_y_window(cfg);
so_plot_phase_background(ax1, cfg, yWindow, false);
draw_rect(ax1, cfg.sourceRectangle, [0.00 0.35 0.85], 'source region');
draw_rect(ax1, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
shown = false;
for r = 1:numel(headRecord.replays)
    pts = split_wrap(so_to_cylinder(headRecord.replays{r}.trajectory));
    if shown
        args = {'HandleVisibility', 'off'};
    else
        args = {'DisplayName', sprintf('noisy replay (\\sigma = %.3g)', headRecord.sigma)};
        shown = true;
    end
    plot(ax1, pts(1, :), pts(2, :), '-', 'Color', [0.96 0.72 0.50], ...
        'LineWidth', 0.6, args{:});
end
detPts = deterministic_polyline(result);
plot(ax1, detPts(1, :), detPts(2, :), '-', 'Color', [0.05 0.05 0.05], ...
    'LineWidth', 1.6, 'DisplayName', 'deterministic schedule');
axis(ax1, [0 1 yWindow]);
xlabel(ax1, 'x  (mod 1)'); ylabel(ax1, 'y');
title(ax1, sprintf('(a) %d realisations at \\sigma = %.3g \\delta', ...
    numel(headRecord.replays), headRecord.factor));
grid(ax1, 'on'); legend(ax1, 'Location', 'best'); hold(ax1, 'off');

% ---- (b) endpoints -------------------------------------------------------
ax2 = nexttile(layout, 2);
hold(ax2, 'on');
draw_rect(ax2, cfg.targetRectangle, [0.85 0.10 0.10], 'target region');
ends = cell2mat(cellfun(@(r) so_to_cylinder(r.finalState), headRecord.replays, ...
    'UniformOutput', false));
inside = cellfun(@(r) r.targetContained, headRecord.replays);
if any(inside)
    plot(ax2, ends(1, inside), ends(2, inside), 'o', 'Color', [0.00 0.35 0.85], ...
        'MarkerFaceColor', [0.60 0.75 0.95], 'MarkerSize', 5, ...
        'DisplayName', 'endpoint in target');
end
if any(~inside)
    plot(ax2, ends(1, ~inside), ends(2, ~inside), 'x', 'Color', [0.85 0.10 0.10], ...
        'MarkerSize', 8, 'LineWidth', 1.3, 'DisplayName', 'endpoint outside');
end
zf = so_to_cylinder(result.finalState);
plot(ax2, zf(1), zf(2), 'p', 'Color', [0.45 0.00 0.55], ...
    'MarkerFaceColor', [0.45 0.00 0.55], 'MarkerSize', 12, ...
    'DisplayName', 'deterministic endpoint');
pad = 0.7 * max(cfg.targetRectangle.xMax - cfg.targetRectangle.xMin, ...
    cfg.targetRectangle.yMax - cfg.targetRectangle.yMin);
axis(ax2, [cfg.targetRectangle.xMin - pad, cfg.targetRectangle.xMax + pad, ...
    cfg.targetRectangle.yMin - pad, cfg.targetRectangle.yMax + pad]);
xlabel(ax2, 'x  (mod 1)'); ylabel(ax2, 'y');
title(ax2, sprintf('(b) endpoints, %d of %d contained', sum(inside), numel(inside)));
grid(ax2, 'on'); legend(ax2, 'Location', 'best'); hold(ax2, 'off');

% ---- (c) containment fraction -------------------------------------------
ax3 = nexttile(layout, 3);
f = levels.sigmaOverDelta;
plotx = f;
% sigma = 0 cannot be drawn on a log axis; park it a decade below the
% smallest real level so the control point stays visible.
if any(f > 0)
    plotx(plotx == 0) = min(f(f > 0)) / 10;
else
    plotx(plotx == 0) = 1e-6;
end
semilogx(ax3, plotx, 100 * levels.containedFraction, 'o-', ...
    'Color', [0.00 0.35 0.85], 'MarkerFaceColor', [0.00 0.35 0.85], ...
    'LineWidth', 1.4);
hold(ax3, 'on');
xline(ax3, noise.headlineFactor, '--', 'published level', ...
    'Color', [0.85 0.35 0.00], 'LineWidth', 1.4);
ylim(ax3, [-5 105]);
hold(ax3, 'off');
xlabel(ax3, '\sigma / \delta   (leftmost point is \sigma = 0)');
ylabel(ax3, 'realisations reaching target  (%)');
title(ax3, '(c) how much noise the open-loop schedule tolerates');
grid(ax3, 'on');

% ---- (d) endpoint displacement ------------------------------------------
ax4 = nexttile(layout, 4);
med = levels.medianDisplacement;
loglog(ax4, plotx, max(med, eps), 'o-', 'Color', [0.00 0.35 0.85], ...
    'MarkerFaceColor', [0.00 0.35 0.85], 'LineWidth', 1.4, ...
    'DisplayName', 'median');
hold(ax4, 'on');
loglog(ax4, plotx, max(levels.maxDisplacement, eps), '^--', ...
    'Color', [0.85 0.35 0.00], 'DisplayName', 'worst case');
targetScale = min(cfg.targetRectangle.xMax - cfg.targetRectangle.xMin, ...
    cfg.targetRectangle.yMax - cfg.targetRectangle.yMin);
yline(ax4, targetScale, ':', 'target size', 'Color', [0.85 0.10 0.10], ...
    'LineWidth', 1.4);
hold(ax4, 'off');
xlabel(ax4, '\sigma / \delta');
ylabel(ax4, 'endpoint displacement');
title(ax4, '(d) endpoint shift versus noise amplitude');
grid(ax4, 'on'); legend(ax4, 'Location', 'best');

title(layout, sprintf(['Noise robustness, %s case: deterministic schedule replayed ' ...
    'open loop, \\delta = %g'], char(cfg.caseName), noise.delta));

if saveFigure
    exportCfg = cfg;
    exportCfg.figureDirectory = fullfile(cfg.outputDirectory, 'noise', 'figures');
    so_export_figure(fig, exportCfg, 'noise_validation');
end
end

% ======================================================================= %

function pts = deterministic_polyline(result)
pts = zeros(2, 0);
for i = 1:numel(result.executionSegments)
    seg = result.executionSegments(i);
    pts = [pts, [NaN; NaN], so_to_cylinder(seg.path)]; %#ok<AGROW>
end
pts = split_wrap(pts);
end

function draw_rect(ax, rect, color, label)
if rect.xMax >= rect.xMin
    xs = [rect.xMin rect.xMax rect.xMax rect.xMin rect.xMin];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
else
    xs = [rect.xMin 1 1 rect.xMin rect.xMin NaN 0 rect.xMax rect.xMax 0 0];
    ys = [rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin NaN ...
        rect.yMin rect.yMin rect.yMax rect.yMax rect.yMin];
end
plot(ax, xs, ys, '-', 'Color', color, 'LineWidth', 1.9, 'DisplayName', label);
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
