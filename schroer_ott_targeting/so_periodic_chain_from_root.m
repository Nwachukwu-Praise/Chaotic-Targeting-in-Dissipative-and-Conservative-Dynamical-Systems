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

chain.eigen = local_chain_eigen(pointsLifted, p, tr, cfg);
chain.stableEigenvalue = chain.eigen.stableEigenvalue;
chain.unstableEigenvalue = chain.eigen.unstableEigenvalue;
end

function eigenData = local_chain_eigen(pointsLifted, p, tr, cfg)
%LOCAL_CHAIN_EIGEN Monodromy spectrum and per-phase-point eigendirections.
%
%   The monodromy matrix is area preserving, so det M = 1 and the two
%   eigenvalues satisfy lambda_+ * lambda_- = 1 with lambda_+ + lambda_- =
%   trace.  For |trace| > 2 both are real and are obtained in closed form,
%   which is more accurate near the thresholds than calling eig on the
%   accumulated product.  For |trace| < 2 the pair is complex conjugate on
%   the unit circle and no stable/unstable splitting exists.
disc = tr^2 - 4;
if disc > 0
    root = sqrt(disc);
    lambdaA = 0.5 * (tr + root);
    lambdaB = 0.5 * (tr - root);
    if abs(lambdaA) >= abs(lambdaB)
        unstableEigenvalue = lambdaA;
        stableEigenvalue = lambdaB;
    else
        unstableEigenvalue = lambdaB;
        stableEigenvalue = lambdaA;
    end
    isHyperbolic = true;
else
    lambdaA = 0.5 * (tr + 1i * sqrt(max(-disc, 0)));
    unstableEigenvalue = lambdaA;
    stableEigenvalue = conj(lambdaA);
    isHyperbolic = false;
end

eigenData.stableEigenvalue = stableEigenvalue;
eigenData.unstableEigenvalue = unstableEigenvalue;
eigenData.absStableEigenvalue = abs(stableEigenvalue);
eigenData.absUnstableEigenvalue = abs(unstableEigenvalue);
eigenData.isHyperbolic = isHyperbolic;
eigenData.lyapunovExponentPerIterate = log(abs(unstableEigenvalue)) / p;

% Per-phase-point monodromy and local invariant directions.  Starting the
% product at phase point i gives the linearisation of the return map at
% that point, whose eigenvectors are the local stable/unstable directions.
eigenData.monodromyByPhasePoint = zeros(2, 2, p);
eigenData.stableDirection = nan(2, p);
eigenData.unstableDirection = nan(2, p);
eigenData.directionResidual = nan(1, p);
for i = 1:p
    Mi = so_jacobian_product(pointsLifted(:, i), p, cfg, 1);
    eigenData.monodromyByPhasePoint(:, :, i) = Mi;
    if ~isHyperbolic
        continue;
    end
    vUnstable = local_eigenvector(Mi, unstableEigenvalue);
    vStable = local_eigenvector(Mi, stableEigenvalue);
    eigenData.unstableDirection(:, i) = vUnstable;
    eigenData.stableDirection(:, i) = vStable;
    eigenData.directionResidual(i) = max( ...
        norm(Mi * vUnstable - unstableEigenvalue * vUnstable), ...
        norm(Mi * vStable - stableEigenvalue * vStable));
end
end

function v = local_eigenvector(M, lambda)
%LOCAL_EIGENVECTOR Unit eigenvector of a 2x2 matrix for a known eigenvalue.
%
%   (M - lambda I) is rank one for a simple eigenvalue, so an eigenvector is
%   obtained by rotating whichever row has the larger norm.  This avoids the
%   ordering ambiguity of eig and keeps the direction tied to the eigenvalue
%   that was selected from the trace.
A = M - lambda * eye(2);
r1 = A(1, :);
r2 = A(2, :);
if norm(r1) >= norm(r2)
    v = [-r1(2); r1(1)];
else
    v = [-r2(2); r2(1)];
end
n = norm(v);
if n == 0 || ~isfinite(n)
    v = [NaN; NaN];
    return;
end
v = v / n;
if v(2) < 0 || (v(2) == 0 && v(1) < 0)
    v = -v;   % declared sign convention: non-negative y component
end
end
