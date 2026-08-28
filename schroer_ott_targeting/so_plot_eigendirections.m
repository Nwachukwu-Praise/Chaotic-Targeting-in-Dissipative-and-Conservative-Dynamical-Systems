function fig = so_plot_eigendirections(chains, cfg, arrowScale)
%SO_PLOT_EIGENDIRECTIONS Local stable/unstable eigendirections of UPO chains.
%
%   fig = SO_PLOT_EIGENDIRECTIONS(chains, cfg) draws, at every phase point
%   of every direct-hyperbolic chain supplied, the two invariant directions
%   of the linearised return map DF^p at that point.
%
%   These are the linear approximations to the stable and unstable
%   manifolds.  They are shown here as a diagnostic of local geometry: they
%   indicate how each resonance organises transport, which is why the route
%   is built from direct-hyperbolic chains only.  They are not used as an
%   explicit pass network and no switch is justified by proximity to them.
%
%   x is displayed modulo one; y is the unwrapped transport coordinate.

if nargin < 3 || isempty(arrowScale)
    arrowScale = 0.045;
end

hyperbolic = chains(strcmp({chains.classification}, 'direct-hyperbolic'));
if isempty(hyperbolic)
    error('SchroerOtt:NoHyperbolicChains', ...
        'No direct-hyperbolic chain was supplied to so_plot_eigendirections.');
end

unstableColour = [0.85 0.15 0.15];
stableColour = [0.15 0.30 0.80];

fig = figure('Name', 'UPO eigendirections', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 17 10]);
layout = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
hold(ax, 'on');

xU = []; yU = []; uU = []; vU = [];
xS = []; yS = []; uS = []; vS = [];
for i = 1:numel(hyperbolic)
    chain = hyperbolic(i);
    eigenData = local_eigen_data(chain, cfg);
    pts = chain.pointsLifted;
    for j = 1:chain.period
        px = so_wrap_x(pts(1, j));
        py = pts(2, j);
        du = eigenData.unstableDirection(:, j);
        ds = eigenData.stableDirection(:, j);
        if all(isfinite(du))
            xU(end + 1) = px; yU(end + 1) = py; %#ok<AGROW>
            uU(end + 1) = du(1); vU(end + 1) = du(2); %#ok<AGROW>
        end
        if all(isfinite(ds))
            xS(end + 1) = px; yS(end + 1) = py; %#ok<AGROW>
            uS(end + 1) = ds(1); vS(end + 1) = ds(2); %#ok<AGROW>
        end
    end
end

% Double-headed arrows: an eigendirection is a line, not an orientation.
hU = quiver(ax, [xU, xU], [yU, yU], arrowScale * [uU, -uU], arrowScale * [vU, -vU], 0, ...
    'Color', unstableColour, 'LineWidth', 0.9, 'MaxHeadSize', 0.5, ...
    'DisplayName', 'unstable direction');
hS = quiver(ax, [xS, xS], [yS, yS], arrowScale * [uS, -uS], arrowScale * [vS, -vS], 0, ...
    'Color', stableColour, 'LineWidth', 0.9, 'MaxHeadSize', 0.5, ...
    'DisplayName', 'stable direction');
plot(ax, xU, yU, 'k.', 'MarkerSize', 6, 'HandleVisibility', 'off');

handles = [hU, hS];
if isfield(cfg, 'sourceRectangle')
    handles(end + 1) = local_rectangle(ax, cfg.sourceRectangle, [0.10 0.35 0.90], 'source');
end
if isfield(cfg, 'targetRectangle')
    handles(end + 1) = local_rectangle(ax, cfg.targetRectangle, [0.90 0.15 0.15], 'target');
end

yAll = [yU, yS];
if isfield(cfg, 'sourceRectangle')
    yAll = [yAll, cfg.sourceRectangle.yMin, cfg.sourceRectangle.yMax];
end
if isfield(cfg, 'targetRectangle')
    yAll = [yAll, cfg.targetRectangle.yMin, cfg.targetRectangle.yMax];
end
pad = 0.06 * max(max(yAll) - min(yAll), 0.1);
axis(ax, [0 1 min(yAll) - pad, max(yAll) + pad]);

xlabel(ax, 'x mod 1');
ylabel(ax, 'y (unwrapped)');
title(ax, 'Stable and unstable eigendirections of direct-hyperbolic chains');
grid(ax, 'on');
box(ax, 'on');
legend(ax, handles, 'Location', 'best');
hold(ax, 'off');

if isfield(cfg, 'saveFigures') && cfg.saveFigures && isfield(cfg, 'figureDirectory')
    if ~exist(cfg.figureDirectory, 'dir')
        mkdir(cfg.figureDirectory);
    end
    base = fullfile(cfg.figureDirectory, 'upo_eigendirections');
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
    exportgraphics(fig, [base '.eps'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    saveas(fig, [base '.svg']);
end
end

% -------------------------------------------------------------------------

function eigenData = local_eigen_data(chain, cfg)
%LOCAL_EIGEN_DATA Use cached eigen data, or recompute for an old catalogue.
if isfield(chain, 'eigen') && ~isempty(chain.eigen) && ...
        isfield(chain.eigen, 'unstableDirection')
    eigenData = chain.eigen;
    return;
end
p = chain.period;
tr = chain.trace;
root = sqrt(max(tr^2 - 4, 0));
a = 0.5 * (tr + root);
b = 0.5 * (tr - root);
if abs(a) >= abs(b)
    lambdaU = a; lambdaS = b;
else
    lambdaU = b; lambdaS = a;
end
eigenData.unstableDirection = nan(2, p);
eigenData.stableDirection = nan(2, p);
for j = 1:p
    Mi = so_jacobian_product(chain.pointsLifted(:, j), p, cfg, 1);
    eigenData.unstableDirection(:, j) = local_eigenvector(Mi, lambdaU);
    eigenData.stableDirection(:, j) = local_eigenvector(Mi, lambdaS);
end
end

function v = local_eigenvector(M, lambda)
A = M - lambda * eye(2);
r1 = A(1, :);
r2 = A(2, :);
if norm(r1) >= norm(r2)
    v = [-r1(2); r1(1)];
else
    v = [-r2(2); r2(1)];
end
n = norm(v);
if n == 0 || ~isfinite(n)
    v = [NaN; NaN];
    return;
end
v = v / n;
if v(2) < 0 || (v(2) == 0 && v(1) < 0)
    v = -v;
end
end

function h = local_rectangle(ax, rect, colour, label)
x = [rect.xMin, rect.xMax, rect.xMax, rect.xMin, rect.xMin];
y = [rect.yMin, rect.yMin, rect.yMax, rect.yMax, rect.yMin];
h = plot(ax, x, y, '-', 'Color', colour, 'LineWidth', 1.4, 'DisplayName', label);
end
