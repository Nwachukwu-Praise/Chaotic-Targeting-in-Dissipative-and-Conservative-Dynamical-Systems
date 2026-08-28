function z = so_rectangle_center(rect)
%SO_RECTANGLE_CENTER Center of a rectangle in its declared x/y lift.
if rect.xMax >= rect.xMin
    x = 0.5 * (rect.xMin + rect.xMax);
else
    x = so_wrap_x(0.5 * (rect.xMin + rect.xMax + 1));
end
y = 0.5 * (rect.yMin + rect.yMax);
z = [x; y];
end

