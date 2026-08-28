function [root, residualNorm, converged] = so_newton_periodic_orbit(seed, p, m, cfg)
%SO_NEWTON_PERIODIC_ORBIT Newton solve F_lift^p(z)-z-[m;0]=0.
z = seed(:);
converged = false;
residualNorm = Inf;
for iter = 1:cfg.orbit.newtonMaxIterations
    [r, A] = residual_and_jacobian(z, p, m, cfg);
    residualNorm = norm(r);
    if residualNorm < cfg.orbit.newtonTolerance
        converged = true;
        break;
    end
    if rcond(A) < 1e-13
        break;
    end
    step = -A \ r;
    damp = 1;
    accepted = false;
    for ls = 1:14
        zTry = z + damp * step;
        [rTry, ~] = residual_and_jacobian(zTry, p, m, cfg);
        if norm(rTry) < residualNorm
            z = zTry;
            accepted = true;
            break;
        end
        damp = 0.5 * damp;
    end
    if ~accepted
        break;
    end
end
[r, ~] = residual_and_jacobian(z, p, m, cfg);
residualNorm = norm(r);
converged = converged || residualNorm < cfg.orbit.residualTolerance;
root = z;
root(1) = so_wrap_x(root(1));
end

function [r, A] = residual_and_jacobian(z, p, m, cfg)
zp = so_iterate(z, p, cfg, 1, false);
M = so_jacobian_product(z, p, cfg, 1);
r = zp - z - [m; 0];
A = M - eye(2);
end

