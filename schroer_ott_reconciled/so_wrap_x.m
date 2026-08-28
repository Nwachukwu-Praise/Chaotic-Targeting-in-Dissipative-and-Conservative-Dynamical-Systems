function x = so_wrap_x(x)
%SO_WRAP_X Wrap x to [0,1).
x = x - floor(x);
x(x >= 1) = x(x >= 1) - 1;
end

