function hits = so_find_polyline_intersections_exhaustive(curveA, curveB, cfg)
%SO_FIND_POLYLINE_INTERSECTIONS_EXHAUSTIVE Reference checker for tests.
segA = so_curve_segments(curveA);
segB = so_curve_segments(curveB);
nA = size(segA.a, 2);
nB = size(segB.a, 2);
hits = struct('segmentA', {}, 'segmentB', {}, 'shiftB', {}, 'alpha', {}, ...
    'beta', {}, 'point', {}, 'sA', {}, 'sB', {});
if nA == 0 || nB == 0
    return;
end
minAX = min(min(segA.a(1, :), segA.b(1, :)));
maxAX = max(max(segA.a(1, :), segA.b(1, :)));
minBX = min(min(segB.a(1, :), segB.b(1, :)));
maxBX = max(max(segB.a(1, :), segB.b(1, :)));
shifts = floor(minAX - maxBX) - 1 : ceil(maxAX - minBX) + 1;
count = 0;
for i = 1:nA
    for j = 1:nB
        for sh = shifts
            h = so_segment_intersection(segA.a(:, i), segA.b(:, i), ...
                segB.a(:, j) + [sh; 0], segB.b(:, j) + [sh; 0], cfg.intersectionTolerance);
            if h.success
                count = count + 1;
                hits(count).segmentA = i; %#ok<AGROW>
                hits(count).segmentB = j; %#ok<AGROW>
                hits(count).shiftB = sh; %#ok<AGROW>
                hits(count).alpha = h.alpha; %#ok<AGROW>
                hits(count).beta = h.beta; %#ok<AGROW>
                hits(count).point = h.point; %#ok<AGROW>
                hits(count).sA = segA.sA(i) + h.alpha * (segA.sB(i) - segA.sA(i)); %#ok<AGROW>
                hits(count).sB = segB.sA(j) + h.beta * (segB.sB(j) - segB.sA(j)); %#ok<AGROW>
            end
        end
    end
end
end

