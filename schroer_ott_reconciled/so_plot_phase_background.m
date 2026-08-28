function h = so_plot_phase_background(ax, cfg, yWindow, showInLegend)
%SO_PLOT_PHASE_BACKGROUND Draw the uncontrolled phase portrait behind a plot.
%
% Uses line/marker rather than scatter: for the ~500k points a readable
% portrait needs, scatter builds one graphics primitive per point and is
% unusably slow, while a single line object with '.' markers is immediate.
if nargin < 4
    showInLegend = true;
end
h = gobjects(0);
if ~cfg.background.enable
    return;
end
points = so_phase_portrait(cfg, yWindow);
if isempty(points)
    return;
end
h = line(ax, points(1, :), points(2, :), ...
    'LineStyle', 'none', 'Marker', '.', ...
    'MarkerSize', cfg.background.markerSize, ...
    'Color', cfg.background.color);
if showInLegend
    set(h, 'DisplayName', 'uncontrolled phase portrait');
else
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
end
uistack(h, 'bottom');
end
