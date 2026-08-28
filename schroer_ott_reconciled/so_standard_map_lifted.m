function zNext = so_standard_map_lifted(z, cfg)
%SO_STANDARD_MAP_LIFTED Lifted normalized standard map.
%
% x is allowed to wind; y is the unbounded transport coordinate.
k = cfg.k;
yNext = z(2, :) - k / (2 * pi) .* sin(2 * pi .* z(1, :));
xNext = z(1, :) + yNext;
zNext = [xNext; yNext];
end

