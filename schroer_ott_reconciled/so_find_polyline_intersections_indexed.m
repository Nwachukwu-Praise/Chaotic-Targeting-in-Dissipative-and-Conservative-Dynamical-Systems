function [hits, stats] = so_find_polyline_intersections_indexed(curveA, curveB, cfg)
%SO_FIND_POLYLINE_INTERSECTIONS_INDEXED Uniform-grid segment intersection.
%
% x-periodicity is handled by adding integer-shifted copies of curveB
% segments that can overlap the lifted x-range of curveA.  y is never
% wrapped.
segA = so_curve_segments(curveA);
segB = so_curve_segments(curveB);
nA = size(segA.a, 2);
nB = size(segB.a, 2);
stats.indexedQueries = 1;
stats.exactTests = 0;
stats.naiveComparisons = nA * nB;
stats.naiveAvoided = 0;

if nA == 0 || nB == 0
    hits = empty_hits();
    stats.naiveAvoided = stats.naiveComparisons;
    return;
end

if exist('polyxpoly', 'file') == 2
    [hits, stats] = fast_polyxpoly_path(segA, segB, nA, nB, curveA, curveB, cfg, stats);
    return;
end

minAX = min(min(segA.a(1, :), segA.b(1, :)));
maxAX = max(max(segA.a(1, :), segA.b(1, :)));
minBX = min(min(segB.a(1, :), segB.b(1, :)));
maxBX = max(max(segB.a(1, :), segB.b(1, :)));
shifts = floor(minAX - maxBX) - 1 : ceil(maxAX - minBX) + 1;

cellSize = max([curveA.hMax, curveB.hMax, 1e-3]);
boxes = [];
copyIndex = 0;
for sh = shifts
    for j = 1:nB
        copyIndex = copyIndex + 1;
        a = segB.a(:, j) + [sh; 0];
        b = segB.b(:, j) + [sh; 0];
        boxes(copyIndex).a = a; %#ok<AGROW>
        boxes(copyIndex).b = b; %#ok<AGROW>
        boxes(copyIndex).segmentIndex = j; %#ok<AGROW>
        boxes(copyIndex).shift = sh; %#ok<AGROW>
        boxes(copyIndex).xmin = min(a(1), b(1)); %#ok<AGROW>
        boxes(copyIndex).xmax = max(a(1), b(1)); %#ok<AGROW>
        boxes(copyIndex).ymin = min(a(2), b(2)); %#ok<AGROW>
        boxes(copyIndex).ymax = max(a(2), b(2)); %#ok<AGROW>
    end
end

grid = containers.Map('KeyType', 'char', 'ValueType', 'any');
for j = 1:numel(boxes)
    ix = floor(boxes(j).xmin / cellSize) : floor(boxes(j).xmax / cellSize);
    iy = floor(boxes(j).ymin / cellSize) : floor(boxes(j).ymax / cellSize);
    for a = ix
        for b = iy
            key = sprintf('%d,%d', a, b);
            if isKey(grid, key)
                grid(key) = [grid(key), j];
            else
                grid(key) = j;
            end
        end
    end
end

hitList = empty_hits();
hitCount = 0;
tested = false(nA, numel(boxes));
for i = 1:nA
    a = segA.a(:, i);
    b = segA.b(:, i);
    xmin = min(a(1), b(1));
    xmax = max(a(1), b(1));
    ymin = min(a(2), b(2));
    ymax = max(a(2), b(2));
    ix = floor(xmin / cellSize) : floor(xmax / cellSize);
    iy = floor(ymin / cellSize) : floor(ymax / cellSize);
    candidates = [];
    for gx = ix
        for gy = iy
            key = sprintf('%d,%d', gx, gy);
            if isKey(grid, key)
                candidates = [candidates, grid(key)]; %#ok<AGROW>
            end
        end
    end
    candidates = unique(candidates(:)).';
    for cidx = candidates
        if tested(i, cidx)
            continue;
        end
        tested(i, cidx) = true;
        box = boxes(cidx);
        if xmax < box.xmin - cfg.intersectionTolerance || xmin > box.xmax + cfg.intersectionTolerance || ...
                ymax < box.ymin - cfg.intersectionTolerance || ymin > box.ymax + cfg.intersectionTolerance
            continue;
        end
        stats.exactTests = stats.exactTests + 1;
        h = so_segment_intersection(a, b, box.a, box.b, cfg.intersectionTolerance);
        if h.success
            hitCount = hitCount + 1;
            hitList(hitCount).segmentA = i; %#ok<AGROW>
            hitList(hitCount).segmentB = box.segmentIndex; %#ok<AGROW>
            hitList(hitCount).shiftB = box.shift; %#ok<AGROW>
            hitList(hitCount).alpha = h.alpha; %#ok<AGROW>
            hitList(hitCount).beta = h.beta; %#ok<AGROW>
            hitList(hitCount).point = h.point; %#ok<AGROW>
            hitList(hitCount).sA = segA.sA(i) + h.alpha * (segA.sB(i) - segA.sA(i)); %#ok<AGROW>
            hitList(hitCount).sB = segB.sA(box.segmentIndex) + h.beta * ...
                (segB.sB(box.segmentIndex) - segB.sA(box.segmentIndex)); %#ok<AGROW>
        end
    end
end

hits = hitList;
stats.naiveAvoided = max(0, stats.naiveComparisons - stats.exactTests);
end

function hits = empty_hits()
hits = struct('segmentA', {}, 'segmentB', {}, 'shiftB', {}, 'alpha', {}, ...
    'beta', {}, 'point', {}, 'sA', {}, 'sB', {});
end

function [hits, stats] = fast_polyxpoly_path(segA, segB, nA, nB, curveA, curveB, cfg, stats)
minAX = min(min(segA.a(1, :), segA.b(1, :)));
maxAX = max(max(segA.a(1, :), segA.b(1, :)));
minBX = min(min(segB.a(1, :), segB.b(1, :)));
maxBX = max(max(segB.a(1, :), segB.b(1, :)));
shifts = floor(minAX - maxBX) - 1 : ceil(maxAX - minBX) + 1;
[xA, yA, mapA] = segment_polyline_vectors(segA, 0);
hits = empty_hits();
hitCount = 0;
for sh = shifts
    [xB, yB, mapB] = segment_polyline_vectors(segB, sh);
    [xi, yi, ii] = polyxpoly(xA, yA, xB, yB);
    stats.exactTests = stats.exactTests + size(ii, 1);
    for k = 1:numel(xi)
        ia = mapA(ii(k, 1));
        ib = mapB(ii(k, 2));
        if ia == 0 || ib == 0
            continue;
        end
        alpha = segment_fraction(segA.a(:, ia), segA.b(:, ia), [xi(k); yi(k)]);
        beta = segment_fraction(segB.a(:, ib) + [sh; 0], segB.b(:, ib) + [sh; 0], [xi(k); yi(k)]);
        hitCount = hitCount + 1;
        hits(hitCount).segmentA = ia; %#ok<AGROW>
        hits(hitCount).segmentB = ib; %#ok<AGROW>
        hits(hitCount).shiftB = sh; %#ok<AGROW>
        hits(hitCount).alpha = alpha; %#ok<AGROW>
        hits(hitCount).beta = beta; %#ok<AGROW>
        hits(hitCount).point = [xi(k); yi(k)]; %#ok<AGROW>
        hits(hitCount).sA = segA.sA(ia) + alpha * (segA.sB(ia) - segA.sA(ia)); %#ok<AGROW>
        hits(hitCount).sB = segB.sA(ib) + beta * (segB.sB(ib) - segB.sA(ib)); %#ok<AGROW>
    end
end
stats.naiveAvoided = max(0, stats.naiveComparisons - stats.exactTests);
% Preserve the diagnostic status as an indexed query.  polyxpoly supplies
% segment indices directly, and all accepted hits still pass Newton
% refinement through the original parametrisations.
stats.implementation = "polyxpoly-segment-indexed";
end

function [x, y, segmentMap] = segment_polyline_vectors(seg, shift)
n = size(seg.a, 2);
x = NaN(1, 3 * n);
y = NaN(1, 3 * n);
segmentMap = zeros(1, 3 * n);
for i = 1:n
    pos = 3 * i - 2;
    a = seg.a(:, i) + [shift; 0];
    b = seg.b(:, i) + [shift; 0];
    x(pos:pos + 1) = [a(1), b(1)];
    y(pos:pos + 1) = [a(2), b(2)];
    segmentMap(pos) = i;
end
end

function f = segment_fraction(a, b, p)
v = b - a;
den = dot(v, v);
if den <= eps
    f = 0;
else
    f = max(0, min(1, dot(p - a, v) / den));
end
end
