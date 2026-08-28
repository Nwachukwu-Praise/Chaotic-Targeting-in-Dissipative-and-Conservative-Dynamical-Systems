function dx = so_wrap_diff_x(dx)
%SO_WRAP_DIFF_X Periodic x difference in [-1/2,1/2).
dx = dx - round(dx);
end

