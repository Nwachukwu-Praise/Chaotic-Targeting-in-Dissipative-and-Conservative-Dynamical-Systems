function chain = so_periodic_chain_from_root(root, p, m, residual, cfg)
%SO_PERIODIC_CHAIN_FROM_ROOT Calculate orbit, trace, residue, class.
pointsLifted = zeros(2, p);
z = root(:);
for i = 1:p
    pointsLifted(:, i) = z;
    z = so_standard_map_lifted(z, cfg);
end
pointsCylinder = so_to_cylinder(pointsLifted);
pointsCylinder(2, :) = pointsLifted(2, :);
M = so_jacobian_product(root(:), p, cfg, 1);
tr = trace(M);
residue = (2 - tr) / 4;
classification = so_classify_trace(tr, cfg.orbit.traceTolerance);

lower = false;
for q = 1:(p - 1)
    if mod(p, q) ~= 0
        continue;
    end
    mq = m * q / p;
    if abs(mq - round(mq)) > 1e-12
        continue;
    end
    rq = so_iterate(root(:), q, cfg, 1, false) - root(:) - [round(mq); 0];
    if norm(rq) < 1e-8
        lower = true;
        break;
    end
end

chain.id = '';
chain.period = p;
chain.winding = m;
chain.omega = m / p;
chain.root = root(:);
chain.pointsLifted = pointsLifted;
chain.pointsCylinder = pointsCylinder;
chain.monodromy = M;
chain.trace = tr;
chain.residue = residue;
chain.absResidue = abs(residue);
chain.residual = residual;
chain.classification = classification;
chain.isLowerPeriod = lower;
chain.averageY = mean(pointsLifted(2, :));
chain.averageYError = abs(chain.averageY - chain.omega);
chain.samePMCount = NaN;
end
