function ecr_plot_radius(model, level, fname)
%ECR_PLOT_RADIUS  The inter-data distance histogram behind the cluster radius.
%
%   ECR_PLOT_RADIUS(MODEL,LEVEL)        histogram of control region T_LEVEL
%   ECR_PLOT_RADIUS(MODEL,LEVEL,FNAME)  also save it as PNG
%
%   Section 3.3 of the paper picks the clustering radius as "the first
%   minimum of the histogram of the inter-data distances by visual analysis".
%   ECR_CHOOSE_RADIUS automates that choice; this figure lets you check it by
%   eye, which is what the paper actually prescribes.

if nargin < 2 || isempty(level), level = 1; end
k = find([model.levels.index] == level);
if isempty(k)
    error('ecr_plot_radius:noLevel', 'the model has no region T%d', level);
end
L = model.levels(k);
if ~isfield(L,'radiusInfo') || ~isfield(L.radiusInfo,'centers') || ...
        isempty(L.radiusInfo.centers)
    error('ecr_plot_radius:noInfo', ...
          'region T%d was not clustered (ECR-II model?)', level);
end

figure('visible','off');
bar(L.radiusInfo.centers, L.radiusInfo.counts, 'FaceColor', [0.75 0.78 0.85], ...
    'EdgeColor', 'none');
hold on
plot(L.radiusInfo.centers, L.radiusInfo.smooth, 'b-', 'LineWidth', 1.4);
yl = ylim;
plot([L.radius L.radius], yl, 'r--', 'LineWidth', 1.6);
xlabel('inter-data distance');
ylabel('count');
ttl = sprintf('T_%d: chosen radius r = %.4g', level, L.radius);
if isfield(L.radiusInfo,'fallback') && L.radiusInfo.fallback
    ttl = [ttl, '  (fallback: no interior minimum)'];
end
title(ttl);
grid on

if nargin >= 3 && ~isempty(fname)
    print(gcf, fname, '-dpng', '-r150');
end
end
