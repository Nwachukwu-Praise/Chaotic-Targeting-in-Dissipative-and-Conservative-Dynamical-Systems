function M = so_jacobian_product(z0, n, cfg, direction)
%SO_JACOBIAN_PRODUCT Product derivative for n forward/backward iterates.
if nargin < 4
    direction = 1;
end
M = eye(2);
z = z0;
for i = 1:n
    if direction >= 0
        J = so_jacobian(z, cfg);
        M = J * M;
        z = so_standard_map_lifted(z, cfg);
    else
        zPrev = so_standard_map_inverse_lifted(z, cfg);
        Jprev = so_jacobian(zPrev, cfg);
        M = (Jprev \ M);
        z = zPrev;
    end
end
end

