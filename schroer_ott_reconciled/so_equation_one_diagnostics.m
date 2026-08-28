function diagnostics = so_equation_one_diagnostics(stagePlans, cfg)
%SO_EQUATION_ONE_DIAGNOSTICS Order-estimate diagnostics for accepted splits.
rows = {};
delta = 2 * cfg.controlAmplitude;
for i = 1:numel(stagePlans)
    conn = stagePlans(i).provisionalConnection;
    if ~conn.success || isinf(conn.totalIterations)
        continue;
    end
    if strcmp(stagePlans(i).targetKind, "final")
        epsT = min(cfg.targetRectangle.xMax - cfg.targetRectangle.xMin, ...
            cfg.targetRectangle.yMax - cfg.targetRectangle.yMin);
    else
        epsT = 2 * cfg.proxyTargetRadius;
    end
    lambda1 = NaN;
    lambda2 = NaN;
    pred = 0;
    if conn.nForward > 0
        M = so_jacobian_product(conn.controlledInitialState, conn.nForward, cfg, 1);
        lambda1 = (1 / conn.nForward) * log(norm(M * [0; 1]));
        if lambda1 > 0
            pred = pred + (1 / lambda1) * log(cfg.equationOne.L / delta);
        end
    end
    if conn.nBackward > 0
        if strcmp(stagePlans(i).targetKind, "final")
            [targetPoint, tangent, cornerAmbiguous] = rectangle_tangent(conn.targetBoundaryPoint, ...
                cfg.targetRectangle, cfg);
        else
            [targetPoint, tangent, cornerAmbiguous] = circle_tangent(conn.targetBoundaryParameter, ...
                conn.targetBoundaryPoint);
        end
        Mback = so_jacobian_product(targetPoint, conn.nBackward, cfg, -1);
        lambda2 = abs((1 / conn.nBackward) * log(norm(Mback * tangent)));
        if lambda2 > 0
            pred = pred + (1 / lambda2) * log(cfg.equationOne.L / epsT);
        end
    else
        cornerAmbiguous = false;
    end
    ratio = conn.totalIterations / max(pred, eps);
    rows(end + 1, :) = {stagePlans(i).stage, conn.nForward, conn.nBackward, ...
        conn.totalIterations, pred, delta, epsT, cfg.equationOne.L, ...
        lambda1, lambda2, string(conn.targetComponent), ratio, cornerAmbiguous}; %#ok<AGROW>
end
names = {'stage','nForward','nBackward','realisedTau','predictedTau','delta', ...
    'epsilon_t','L','lambda1','absLambda2','targetComponent','discrepancyRatio'};
names{end + 1} = 'cornerAmbiguous';
if isempty(rows)
    diagnostics = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    diagnostics = cell2table(rows, 'VariableNames', names);
end
end

function [point, tangent, cornerAmbiguous] = circle_tangent(s, boundaryPoint)
theta = 2 * pi * s;
tangent = [-sin(theta); cos(theta)];
tangent = tangent / norm(tangent);
point = boundaryPoint(:);
cornerAmbiguous = false;
end

function [point, tangent, cornerAmbiguous] = rectangle_tangent(point, rect, cfg)
rect.yMin = rect.yMin + cfg.transport.targetLiftShift;
rect.yMax = rect.yMax + cfg.transport.targetLiftShift;
zc = so_to_cylinder(point(:));
tol = 1e-7;
onLeft = abs(so_wrap_diff_x(zc(1) - rect.xMin)) < tol;
onRight = abs(so_wrap_diff_x(zc(1) - rect.xMax)) < tol;
onBottom = abs(zc(2) - rect.yMin) < tol;
onTop = abs(zc(2) - rect.yMax) < tol;
cornerAmbiguous = (onLeft || onRight) && (onBottom || onTop);
if onBottom || onTop
    tangent = [1; 0];
else
    tangent = [0; 1];
end
point = point(:);
end
