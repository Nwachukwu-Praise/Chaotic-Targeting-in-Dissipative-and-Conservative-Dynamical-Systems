function refined = so_refine_intersection(sourceComponent, targetComponent, ...
    nForward, nBackward, u0, s0, targetShift, cfg, uBracket, sBracket)
%SO_REFINE_INTERSECTION Bracket-preserving refinement in original parameters.
%
% The polyline search identifies one pair of intersecting curve segments.
% This routine refines only inside those two parameter brackets, so a
% successful Newton step cannot jump to a different crossing.
if nargin < 9 || isempty(uBracket)
    uBracket = sourceComponent.parameterRange;
end
if nargin < 10 || isempty(sBracket)
    sBracket = targetComponent.parameterRange;
end

uRange = sourceComponent.parameterRange;
uLo = max(min(uBracket), min(uRange));
uHi = min(max(uBracket), max(uRange));
if uLo > uHi
    uLo = min(uRange);
    uHi = max(uRange);
end

[sLo, sHi, sStart] = normalized_s_bracket(targetComponent, sBracket, s0);
u = clamp(u0, uLo, uHi);
s = clamp(sStart, sLo, sHi);

best.u = u;
best.s = s;
[best.r, best.zF, best.zB] = residual_at(u, s);
best.norm = norm(best.r);
acceptedNewton = false;

for iter = 1:25
    [r, zF, zB] = residual_at(u, s);
    residual = norm(r);
    if residual < best.norm
        best = make_best(u, s, r, zF, zB);
    end
    if residual <= cfg.intersectionTolerance
        acceptedNewton = true;
        break;
    end

    duCol = finite_difference_column(u, s, 1);
    dsCol = finite_difference_column(u, s, 2);
    A = [duCol, dsCol];
    if any(~isfinite(A(:))) || rcond(A) < 1e-12
        break;
    end

    step = -A \ r;
    damp = 1;
    acceptedStep = false;
    for ls = 1:14
        uTry = clamp(u + damp * step(1), uLo, uHi);
        sTry = clamp(s + damp * step(2), sLo, sHi);
        rTry = residual_at(uTry, sTry);
        if norm(rTry) < residual
            u = uTry;
            s = sTry;
            acceptedStep = true;
            break;
        end
        damp = 0.5 * damp;
    end
    if ~acceptedStep
        break;
    end
end

if best.norm > cfg.intersectionTolerance
    best = local_subdivision_best(best);
end

u = best.u;
s = best.s;
zF = best.zF;
zB = best.zB;
residual = best.norm;
zControlled = sourceComponent.gamma0(u);
plannedPath = so_iterate(zControlled, nForward + nBackward, cfg, 1, true);
finalState = plannedPath(:, end);
targetContained = targetComponent.contains(finalState, cfg.containmentTolerance);

refined.success = (acceptedNewton || residual <= 10 * cfg.intersectionTolerance) && targetContained;
refined.control = u;
refined.boundaryParameter = mod(s, 1);
refined.intersectionPoint = zF;
refined.targetBoundaryPoint = zB - [targetShift; 0];
refined.intersectionResidual = residual;
refined.controlledInitialState = zControlled;
refined.preControlState = sourceComponent.center;
refined.postControlState = zControlled;
refined.plannedPath = plannedPath;
refined.finalState = finalState;
refined.targetContained = targetContained;
refined.bracket = [uLo, uHi, sLo, sHi];
refined.bracketPreserved = u >= uLo - eps(uLo) && u <= uHi + eps(uHi) && ...
    s >= sLo - eps(sLo) && s <= sHi + eps(sHi);

    function [r, zF, zB] = residual_at(uu, ss)
        zF = so_curve_point(sourceComponent, uu, 1, nForward, cfg);
        zB = so_curve_point(targetComponent, ss, -1, nBackward, cfg) + [targetShift; 0];
        zF(1) = so_lift_x_near(so_wrap_x(zF(1)), zB(1));
        r = zF - zB;
    end

    function col = finite_difference_column(uu, ss, whichParameter)
        if whichParameter == 1
            lo = uLo;
            hi = uHi;
            width = hi - lo;
            step = min(max(1e-8, 1e-4 * max(width, cfg.controlAmplitude)), 0.25 * max(width, eps));
            [pMinus, pPlus] = bounded_stencil(uu, lo, hi, step);
            if pPlus == pMinus
                col = [NaN; NaN];
                return;
            end
            rMinus = residual_at(pMinus, ss);
            rPlus = residual_at(pPlus, ss);
        else
            lo = sLo;
            hi = sHi;
            width = hi - lo;
            step = min(max(1e-8, 1e-4 * max(width, 1e-3)), 0.25 * max(width, eps));
            [pMinus, pPlus] = bounded_stencil(ss, lo, hi, step);
            if pPlus == pMinus
                col = [NaN; NaN];
                return;
            end
            rMinus = residual_at(uu, pMinus);
            rPlus = residual_at(uu, pPlus);
        end
        col = (rPlus - rMinus) ./ (pPlus - pMinus);
    end

    function bestOut = local_subdivision_best(bestIn)
        bestOut = bestIn;
        uWindow = [uLo, uHi];
        sWindow = [sLo, sHi];
        for pass = 1:4
            uGrid = linspace(uWindow(1), uWindow(2), 9);
            sGrid = linspace(sWindow(1), sWindow(2), 9);
            for ii = 1:numel(uGrid)
                for jj = 1:numel(sGrid)
                    [rr, zzF, zzB] = residual_at(uGrid(ii), sGrid(jj));
                    nr = norm(rr);
                    if nr < bestOut.norm
                        bestOut = make_best(uGrid(ii), sGrid(jj), rr, zzF, zzB);
                    end
                end
            end
            uHalf = max((uWindow(2) - uWindow(1)) / 8, eps);
            sHalf = max((sWindow(2) - sWindow(1)) / 8, eps);
            uWindow = [clamp(bestOut.u - uHalf, uLo, uHi), clamp(bestOut.u + uHalf, uLo, uHi)];
            sWindow = [clamp(bestOut.s - sHalf, sLo, sHi), clamp(bestOut.s + sHalf, sLo, sHi)];
        end
    end
end

function best = make_best(u, s, r, zF, zB)
best.u = u;
best.s = s;
best.r = r;
best.zF = zF;
best.zB = zB;
best.norm = norm(r);
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function [lo, hi, sStart] = normalized_s_bracket(component, bracket, s0)
if component.isClosed
    rawLo = min(bracket);
    rawHi = max(bracket);
    width = rawHi - rawLo;
    if width <= 0
        rawLo = floor(s0);
        rawHi = rawLo + 1;
    end
    candidates = s0 + (-2:2);
    inside = candidates >= rawLo - 10 * eps(rawLo) & candidates <= rawHi + 10 * eps(rawHi);
    if any(inside)
        sStart = candidates(find(inside, 1, 'first'));
        lo = rawLo;
        hi = rawHi;
    else
        sStart = s0;
        lo = rawLo;
        hi = rawHi;
    end
else
    lo = max(min(bracket), min(component.parameterRange));
    hi = min(max(bracket), max(component.parameterRange));
    sStart = s0;
end
end

function [pMinus, pPlus] = bounded_stencil(p, lo, hi, h)
if p - h >= lo && p + h <= hi
    pMinus = p - h;
    pPlus = p + h;
elseif p + h <= hi
    pMinus = p;
    pPlus = p + h;
elseif p - h >= lo
    pMinus = p - h;
    pPlus = p;
else
    pMinus = lo;
    pPlus = hi;
end
end
