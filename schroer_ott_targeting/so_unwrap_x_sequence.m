function points = so_unwrap_x_sequence(points)
%SO_UNWRAP_X_SEQUENCE Sequentially unwrap the x row of a 2-by-N curve.
if isempty(points)
    return;
end
for j = 2:size(points, 2)
    points(1, j) = so_lift_x_near(points(1, j), points(1, j - 1));
end
end

