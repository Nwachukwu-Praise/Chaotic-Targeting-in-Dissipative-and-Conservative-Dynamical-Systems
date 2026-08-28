function manifolds = so_compute_diagnostic_manifolds(route, cfg)
%SO_COMPUTE_DIAGNOSTIC_MANIFOLDS Non-operational manifold diagnostics.
%
% These branches are used only for a posteriori visualisation of the
% selected switching points.  They do not define the route, do not define
% target acceptance, and do not define switching.
branchLength = 0.006;
branchIterates = 10;
manifolds = repmat(empty_record(), 0, 1);
for c = 1:numel(route.chains)
    chain = route.chains(c);
    for phase = 1:chain.period
        z = chain.pointsLifted(:, phase);
        M = so_jacobian_product(z, chain.period, cfg, 1);
        [V, D] = eig(M);
        eigVals = diag(D);
        [~, unstableIdx] = max(abs(eigVals));
        [~, stableIdx] = min(abs(eigVals));
        vu = real(V(:, unstableIdx));
        vs = real(V(:, stableIdx));
        vu = vu / norm(vu);
        vs = vs / norm(vs);
        for signValue = [-1, 1]
            uComponent = make_branch_component(z, signValue * vu, branchLength, ...
                sprintf('%s-phase-%02d-unstable-%+d', chain.id, phase, signValue));
            sComponent = make_branch_component(z, signValue * vs, branchLength, ...
                sprintf('%s-phase-%02d-stable-%+d', chain.id, phase, signValue));
            uCurve = so_build_curve(uComponent, 1, branchIterates, cfg);
            sCurve = so_build_curve(sComponent, -1, branchIterates, cfg);
            manifolds(end + 1) = make_record(chain, phase, "unstable", signValue, uCurve); %#ok<AGROW>
            manifolds(end + 1) = make_record(chain, phase, "stable", signValue, sCurve); %#ok<AGROW>
        end
    end
end
end

function component = make_branch_component(z, directionVector, lengthScale, id)
component.type = 'manifold-branch-seed';
component.id = id;
component.isClosed = false;
component.parameterRange = [0, 1];
component.characteristicSize = lengthScale;
component.center = z(:);
component.phasePointIndex = NaN;
component.chainID = "";
component.gamma0 = @(s) z(:) + directionVector(:) .* (lengthScale .* reshape(s, 1, []));
component.contains = @(zq, tol) false(1, size(zq, 2)); %#ok<INUSD>
end

function rec = make_record(chain, phase, branchType, signValue, curve)
rec = empty_record();
rec.chainID = string(chain.id);
rec.omega = chain.omega;
rec.phasePointIndex = phase;
rec.branchType = branchType;
rec.sign = signValue;
rec.curve = curve;
rec.resolutionStatus = string(curve.resolutionStatus);
rec.pointCount = curve.pointCount;
rec.maximumGap = curve.maximumGap;
rec.maximumMidpointDeviation = curve.maximumMidpointDeviation;
end

function rec = empty_record()
rec.chainID = "";
rec.omega = NaN;
rec.phasePointIndex = NaN;
rec.branchType = "";
rec.sign = NaN;
rec.curve = [];
rec.resolutionStatus = "";
rec.pointCount = NaN;
rec.maximumGap = NaN;
rec.maximumMidpointDeviation = NaN;
end

