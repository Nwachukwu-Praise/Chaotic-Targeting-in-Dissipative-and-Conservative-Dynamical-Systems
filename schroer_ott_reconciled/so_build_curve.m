function curve = so_build_curve(component, direction, iterateCount, cfg)
%SO_BUILD_CURVE Adaptive fold-aware curve image/preimage.
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

for depth = 1:cfg.curve.maxSubdivisionDepth
    [needs, gapNow, midNow] = intervals_needing_refinement(component, parameters, ...
        direction, iterateCount, cfg, hMax, hMid);
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
    [needs, gapNow, midNow] = intervals_needing_refinement(component, parameters, ...
        direction, iterateCount, cfg, hMax, hMid);
    maxGap = max(maxGap, gapNow);
    maxMid = max(maxMid, midNow);
    if any(needs) && status == "resolved"
        status = "unresolved_subdivision_depth";
    end
end

pointsLifted = so_curve_point(component, parameters, direction, iterateCount, cfg);
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

function [needs, maxGap, maxMid] = intervals_needing_refinement(component, parameters, direction, iterateCount, cfg, hMax, hMid)
n = numel(parameters);
if component.isClosed
    intervalCount = n;
else
    intervalCount = n - 1;
end
needs = false(1, intervalCount);
maxGap = 0;
maxMid = 0;
for i = 1:intervalCount
    s1 = parameters(i);
    if i == n
        s2 = parameters(1) + 1;
    else
        s2 = parameters(i + 1);
    end
    sm = 0.5 * (s1 + s2);
    p1 = so_curve_point(component, s1, direction, iterateCount, cfg);
    p2 = so_curve_point(component, s2, direction, iterateCount, cfg);
    pm = so_curve_point(component, sm, direction, iterateCount, cfg);
    p1(1) = so_wrap_x(p1(1));
    p2(1) = so_lift_x_near(so_wrap_x(p2(1)), p1(1));
    pm(1) = so_lift_x_near(so_wrap_x(pm(1)), 0.5 * (p1(1) + p2(1)));
    gap = norm(p2 - p1);
    midDev = norm(pm - 0.5 * (p1 + p2));
    maxGap = max(maxGap, gap);
    maxMid = max(maxMid, midDev);
    needs(i) = gap > hMax || midDev > hMid;
end
end

function parameters = insert_midpoints(component, parameters, needs)
n = numel(parameters);
newParams = parameters(:).';
for i = 1:numel(needs)
    if ~needs(i)
        continue;
    end
    s1 = parameters(i);
    if i == n && component.isClosed
        s2 = parameters(1) + 1;
    else
        s2 = parameters(i + 1);
    end
    sm = 0.5 * (s1 + s2);
    if component.isClosed
        sm = mod(sm, 1);
    end
    newParams(end + 1) = sm; %#ok<AGROW>
end
if component.isClosed
    parameters = unique(mod(newParams, 1));
else
    parameters = unique(newParams);
end
parameters = sort(parameters);
end

