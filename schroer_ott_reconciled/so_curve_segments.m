function segments = so_curve_segments(curve)
%SO_CURVE_SEGMENTS Return ordered line segments with parameter brackets.
p = curve.pointsLifted;
s = curve.parameters;
n = size(p, 2);
if curve.isClosed
    segmentCount = n;
else
    segmentCount = max(0, n - 1);
end
segments.a = zeros(2, segmentCount);
segments.b = zeros(2, segmentCount);
segments.sA = zeros(1, segmentCount);
segments.sB = zeros(1, segmentCount);
for i = 1:segmentCount
    j = i + 1;
    if j > n
        j = 1;
    end
    a = p(:, i);
    b = p(:, j);
    if curve.isClosed && j == 1
        b(1) = so_lift_x_near(b(1), a(1));
        sB = s(1) + 1;
    else
        b(1) = so_lift_x_near(b(1), a(1));
        sB = s(j);
    end
    segments.a(:, i) = a;
    segments.b(:, i) = b;
    segments.sA(i) = s(i);
    segments.sB(i) = sB;
end
end

