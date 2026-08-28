function zc = so_to_cylinder(z)
%SO_TO_CYLINDER Return (x mod 1, y) without wrapping y.
zc = z;
zc(1, :) = so_wrap_x(zc(1, :));
end

