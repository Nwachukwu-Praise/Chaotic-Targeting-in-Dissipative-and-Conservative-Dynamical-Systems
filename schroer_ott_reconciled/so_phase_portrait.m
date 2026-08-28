function points = so_phase_portrait(cfg, yWindow)
%SO_PHASE_PORTRAIT Uncontrolled phase portrait of the standard map.
%
% Returns a 2-by-N array of (x mod 1, y) iterates of the UNCONTROLLED map,
% for use as the background of every phase-space figure.  This is the
% "map underneath": islands show up as closed curves, the connected chaotic
% region as stipple, and the resonance chains the route uses then sit
% visibly inside their own islands instead of floating on an empty axis.
%
% Nothing here feeds the targeting calculation.  It is a picture of the
% system, computed from the same so_standard_map_lifted used everywhere else.
%
% The result is memoised for the session, keyed on k, the window and the
% seed counts, because the same portrait is drawn on several figures.

persistent cacheKey cachePoints
if nargin < 2 || isempty(yWindow)
    yWindow = [-0.1, 1.1];
end
yWindow = [min(yWindow), max(yWindow)];

nx = cfg.background.xSeedCount;
ny = cfg.background.ySeedCount;
nIter = cfg.background.iterations;
key = [cfg.k, yWindow, nx, ny, nIter];

if ~isempty(cacheKey) && isequal(cacheKey, key)
    points = cachePoints;
    return;
end

xSeeds = linspace(0, 1, nx + 1);
xSeeds(end) = [];
ySeeds = linspace(yWindow(1), yWindow(2), ny);

Z = zeros(2, nx * ny);
idx = 0;
for ix = 1:nx
    for iy = 1:ny
        idx = idx + 1;
        Z(:, idx) = [xSeeds(ix); ySeeds(iy)];
    end
end

X = zeros(numel(Z(1, :)), nIter);
Y = zeros(numel(Z(1, :)), nIter);
for n = 1:nIter
    Z = so_standard_map_lifted(Z, cfg);
    X(:, n) = so_wrap_x(Z(1, :)).';
    Y(:, n) = Z(2, :).';
end

x = X(:).';
y = Y(:).';
keep = y >= yWindow(1) & y <= yWindow(2);
points = [x(keep); y(keep)];

cacheKey = key;
cachePoints = points;
end
