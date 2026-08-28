function J = so_jacobian(z, cfg)
%SO_JACOBIAN Analytic Jacobian of [x_{n+1}; y_{n+1}].
c = cfg.k .* cos(2 * pi .* z(1));
J = [1 - c, 1; -c, 1];
end

