function trajectory = so_iterate(z0, n, cfg, direction, keepHistory)
%SO_ITERATE Iterate lifted map forward or backward.
if nargin < 4 || isempty(direction)
    direction = 1;
end
if nargin < 5
    keepHistory = false;
end
z = z0;
if keepHistory
    trajectory = zeros(2, n + 1);
    trajectory(:, 1) = z0;
else
    trajectory = z0;
end
for i = 1:n
    if direction >= 0
        z = so_standard_map_lifted(z, cfg);
    else
        z = so_standard_map_inverse_lifted(z, cfg);
    end
    if keepHistory
        trajectory(:, i + 1) = z;
    end
end
if ~keepHistory
    trajectory = z;
end
end

