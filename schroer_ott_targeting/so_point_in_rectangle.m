function inside = so_point_in_rectangle(z, rect, tol)
%SO_POINT_IN_RECTANGLE Rectangle containment, periodic only in x.
if nargin < 3
    tol = 0;
end
zc = so_to_cylinder(z);
if rect.xMax >= rect.xMin
    insideX = zc(1, :) >= rect.xMin - tol & zc(1, :) <= rect.xMax + tol;
else
    insideX = zc(1, :) >= rect.xMin - tol | zc(1, :) <= rect.xMax + tol;
end
insideY = zc(2, :) >= rect.yMin - tol & zc(2, :) <= rect.yMax + tol;
inside = insideX & insideY;
end

