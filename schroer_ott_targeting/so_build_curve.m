function [curve, stats] = so_build_curve(component, direction, iterateCount, cfg)
%SO_BUILD_CURVE Adaptive fold-aware curve image/preimage.
%
%   The refinement mathematics is identical to so_build_curve_reference.
%   Only the evaluation strategy differs:
%
%     (a) all interval endpoints and midpoints of a depth pass are evaluated
%         in a single vectorised so_curve_point call instead of three scalar
%         calls per interval;
%     (b) every evaluated parameter is memoised, so a midpoint computed at
%         depth d is reused when it becomes an inserted node at depth d+1,
%         and an interval that has already converged costs no further map
%         iteration when it is retested.
%
%   The parameter set, the resolution status, the reported maximum gap and
%   midpoint deviation, and the returned points are therefore unchanged.
%   test_so_build_curve_equivalence asserts this against the reference.
%
%   The optional second output reports evaluation counts for profiling.

if component.isClosed
    n0 = cfg.curve.initialBoundarySamples;
    parameters = linspace(0, 1, n0 + 1);
    parameters(end) = [];
else
    n0 = cfg.curve.initialControlSamples;
    parameters = linspace(component.parameterRange(1), component.parameterRange(2), n0);
end

hMax = max(cfg.curve.maxGapFraction * component.characteristicSize, 1e-6);
hMid = max(cfg.curve.midpointToleranceFraction * component.characteristicSize, 1e-7);
status = "resolved";
maxGap = 0;
maxMid = 0;
depth = 0;

cache = new_cache();
stats = struct('curvePointEvaluations', 0, 'mapIterations', 0, ...
    'cacheHits', 0, 'depthPasses', 0);

for depth = 1:cfg.curve.maxSubdivisionDepth
    [needs, gapNow, midNow, cache, stats] = intervals_needing_refinement(component, ...
        parameters, direction, iterateCount, cfg, hMax, hMid, cache, stats);
    maxGap = max(maxGap, gapNow);
    maxMid = max(maxMid, midNow);
    if ~any(needs)
        break;
    end
    parameters = insert_midpoints(component, parameters, needs);
    if numel(parameters) > cfg.curve.maxPoints
        status = "unresolved_point_budget";
        parameters = parameters(1:min(numel(parameters), cfg.curve.maxPoints));
        break;
    end
end

if depth >= cfg.curve.maxSubdivisionDepth
    [needs, gapNow, midNow, cache, stats] = intervals_needing_refinement(component, ...
        parameters, direction, iterateCount, cfg, hMax, hMid, cache, stats);
    maxGap = max(maxGap, gapNow);
    maxMid = max(maxMid, midNow);
    if any(needs) && status == "resolved"
        status = "unresolved_subdivision_depth";
    end
end

[pointsLifted, cache, stats] = cached_points(cache, stats, component, parameters, ...
    direction, iterateCount, cfg);
pointsLifted(1, :) = so_wrap_x(pointsLifted(1, :));
pointsLifted = so_unwrap_x_sequence(pointsLifted);

curve.gamma0 = component.gamma0;
curve.component = component;
curve.parameters = parameters;
curve.isClosed = component.isClosed;
curve.direction = direction;
curve.iterateCount = iterateCount;
curve.pointsLifted = pointsLifted;
curve.pointsCylinder = so_to_cylinder(pointsLifted);
curve.resolutionStatus = status;
curve.maximumGap = maxGap;
curve.maximumMidpointDeviation = maxMid;
curve.subdivisionDepth = depth;
curve.pointCount = numel(parameters);
curve.hMax = hMax;
curve.hMid = hMid;
end

% -------------------------------------------------------------------------

function [needs, maxGap, maxMid, cache, stats] = intervals_needing_refinement(component, ...
    parameters, direction, iterateCount, cfg, hMax, hMid, cache, stats)
n = numel(parameters);
if component.isClosed
    intervalCount = n;
    s1 = parameters;
    s2 = [parameters(2:end), parameters(1) + 1];
else
    intervalCount = n - 1;
    s1 = parameters(1:end - 1);
    s2 = parameters(2:end);
end
needs = false(1, max(intervalCount, 0));
maxGap = 0;
maxMid = 0;
if intervalCount < 1
    return;
end
stats.depthPasses = stats.depthPasses + 1;

sm = 0.5 * (s1 + s2);
[P1, cache, stats] = cached_points(cache, stats, component, s1, direction, iterateCount, cfg);
[P2, cache, stats] = cached_points(cache, stats, component, s2, direction, iterateCount, cfg);
[PM, cache, stats] = cached_points(cache, stats, component, sm, direction, iterateCount, cfg);

P1(1, :) = so_wrap_x(P1(1, :));
P2(1, :) = so_lift_x_near(so_wrap_x(P2(1, :)), P1(1, :));
PM(1, :) = so_lift_x_near(so_wrap_x(PM(1, :)), 0.5 * (P1(1, :) + P2(1, :)));

gap = vecnorm(P2 - P1, 2, 1);
midDev = vecnorm(PM - 0.5 * (P1 + P2), 2, 1);

maxGap = max([0, gap]);
maxMid = max([0, midDev]);
needs = gap > hMax | midDev > hMid;
end

% -------------------------------------------------------------------------

function parameters = insert_midpoints(component, parameters, needs)
if component.isClosed
    s1 = parameters;
    s2 = [parameters(2:end), parameters(1) + 1];
else
    s1 = parameters(1:end - 1);
    s2 = parameters(2:end);
end
sm = 0.5 * (s1(needs) + s2(needs));
if component.isClosed
    sm = mod(sm, 1);
    parameters = unique(mod([parameters(:).', sm(:).'], 1));
else
    parameters = unique([parameters(:).', sm(:).']);
end
parameters = sort(parameters);
end

% -------------------------------------------------------------------------

function cache = new_cache()
cache.s = zeros(1, 0);
cache.p = zeros(2, 0);
cache.count = 0;
end

function [pts, cache, stats] = cached_points(cache, stats, component, s, direction, iterateCount, cfg)
%CACHED_POINTS Evaluate F^(+-n)(gamma0(s)) with exact-key memoisation.
%
%   Keys are the raw double parameter values.  Parameters produced by the
%   refinement are reproduced bit-exactly across depth passes (unique and
%   sort return the input doubles unchanged), so exact equality is the
%   correct lookup rule here and no tolerance is involved.
s = s(:).';
if isempty(s)
    pts = zeros(2, 0);
    return;
end

known = cache.s(1:cache.count);
[found, loc] = ismember(s, known);
if ~all(found)
    missing = unique(s(~found));
    newPoints = so_curve_point(component, missing, direction, iterateCount, cfg);
    stats.curvePointEvaluations = stats.curvePointEvaluations + numel(missing);
    stats.mapIterations = stats.mapIterations + numel(missing) * iterateCount;
    cache = append_cache(cache, missing, newPoints);
    known = cache.s(1:cache.count);
    [found, loc] = ismember(s, known);
end
stats.cacheHits = stats.cacheHits + sum(found);
pts = cache.p(:, loc);
end

function cache = append_cache(cache, s, p)
m = numel(s);
capacity = numel(cache.s);
if cache.count + m > capacity
    newCapacity = max(2 * capacity, cache.count + m);
    cache.s(1, newCapacity) = 0;
    cache.p(2, newCapacity) = 0;
end
idx = cache.count + (1:m);
cache.s(idx) = s;
cache.p(:, idx) = p;
cache.count = cache.count + m;
end
