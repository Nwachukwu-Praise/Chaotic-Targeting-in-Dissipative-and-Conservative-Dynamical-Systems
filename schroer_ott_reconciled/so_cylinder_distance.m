function d = so_cylinder_distance(a, b)
%SO_CYLINDER_DISTANCE Euclidean distance on S^1 x R.
dx = so_wrap_diff_x(a(1, :) - b(1, :));
dy = a(2, :) - b(2, :);
d = hypot(dx, dy);
end

