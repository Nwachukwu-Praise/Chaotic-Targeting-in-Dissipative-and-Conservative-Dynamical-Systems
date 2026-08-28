function refined = so_refine_intersection(sourceComponent, targetComponent, ...
    nForward, nBackward, u0, s0, targetShift, cfg)
%SO_REFINE_INTERSECTION Newton refinement in original parameters.
%
% All trial image points are recomputed from gamma0 and map iterates.
uLo = sourceComponent.parameterRange(1);
uHi = sourceComponent.parameterRange(2);
sBase = floor(s0);
sLo = sBase;
sHi = sBase + 1;
u = min(max(u0, uLo), uHi);
s = s0;
ok = false;
residual = Inf;

for iter = 1:25
    [r, zF, zB] = residual_at(u, s);
    residual = norm(r);
    if residual <= cfg.intersectionTolerance
        ok = true;
        break;
    end
    du = max(1e-8, 1e-5 * cfg.controlAmplitude);
    ds = 1e-6;
    rU = residual_at(min(uHi, u + du), s);
    rS = residual_at(u, s + ds);
    A = [(rU - r) ./ max(min(uHi, u + du) - u, eps), (rS - r) ./ ds];
    if rcond(A) < 1e-12
        break;
    end
    step = -A \ r;
    damp = 1;
    accepted = false;
    for ls = 1:12
        uTry = min(max(u + damp * step(1), uLo), uHi);
        sTry = min(max(s + damp * step(2), sLo), sHi);
        rTry = residual_at(uTry, sTry);
        if norm(rTry) < residual
            u = uTry;
            s = sTry;
            accepted = true;
            break;
        end
        damp = 0.5 * damp;
    end
    if ~accepted
        break;
    end
end

[r, zF, zB] = residual_at(u, s);
residual = norm(r);
ok = ok || residual <= 10 * cfg.intersectionTolerance;
zControlled = sourceComponent.gamma0(u);
plannedPath = so_iterate(zControlled, nForward + nBackward, cfg, 1, true);
finalState = plannedPath(:, end);
targetContained = targetComponent.contains(finalState, cfg.containmentTolerance);

refined.success = ok && targetContained;
refined.control = u;
refined.boundaryParameter = mod(s, 1);
refined.intersectionPoint = zF;
refined.targetBoundaryPoint = zB - [targetShift; 0];
refined.intersectionResidual = residual;
refined.controlledInitialState = zControlled;
refined.plannedPath = plannedPath;
refined.finalState = finalState;
refined.targetContained = targetContained;

    function [r, zF, zB] = residual_at(uu, ss)
        zF = so_curve_point(sourceComponent, uu, 1, nForward, cfg);
        zB = so_curve_point(targetComponent, ss, -1, nBackward, cfg) + [targetShift; 0];
        zF(1) = so_lift_x_near(so_wrap_x(zF(1)), zB(1));
        r = zF - zB;
    end
end
