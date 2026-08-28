function zPrev = so_standard_map_inverse_lifted(z, cfg)
%SO_STANDARD_MAP_INVERSE_LIFTED Inverse in the same chosen lift.
k = cfg.k;
xPrev = z(1, :) - z(2, :);
yPrev = z(2, :) + k / (2 * pi) .* sin(2 * pi .* xPrev);
zPrev = [xPrev; yPrev];
end

