function diagnostics = so_equation_one_diagnostics(stagePlans, cfg)
%SO_EQUATION_ONE_DIAGNOSTICS Order-estimate diagnostics for accepted splits.
rows = {};
delta = 2 * cfg.controlAmplitude;
for i = 1:numel(stagePlans)
    conn = stagePlans(i).provisionalConnection;
    if ~conn.success || isinf(conn.totalIterations)
        continue;
    end
    [epsT, epsAlternative, convention, alternativeConvention] = target_size(stagePlans(i), cfg);
    lambda1 = NaN;
    lambda2 = NaN;
    pred = 0;
    predAlternative = 0;
    if conn.nForward > 0
        M = so_jacobian_product(conn.controlledInitialState, conn.nForward, cfg, 1);
        lambda1 = (1 / conn.nForward) * log(norm(M * [0; 1]));
        if lambda1 > 0
            forwardTerm = (1 / lambda1) * log(cfg.equationOne.L / delta);
            pred = pred + forwardTerm;
            predAlternative = predAlternative + forwardTerm;
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
            predAlternative = predAlternative + (1 / lambda2) * ...
                log(cfg.equationOne.L / epsAlternative);
        end
    else
        cornerAmbiguous = false;
    end
    ratio = conn.totalIterations / max(pred, eps);
    rows(end + 1, :) = {stagePlans(i).stage, conn.nForward, conn.nBackward, ...
        conn.totalIterations, pred, predAlternative, delta, epsT, epsAlternative, ...
        string(convention), string(alternativeConvention), cfg.equationOne.L, ...
        string(cfg.equationOne.LConvention), ...
        lambda1, lambda2, string(conn.targetComponent), ratio, cornerAmbiguous}; %#ok<AGROW>
end
names = {'stage','nForward','nBackward','realisedTau','predictedTau','alternativePredictedTau', ...
    'delta','epsilon_t','alternativeEpsilon_t','targetSizeConvention', ...
    'alternativeTargetSizeConvention','L','LConvention','lambda1','absLambda2', ...
    'targetComponent','discrepancyRatio'};
names{end + 1} = 'cornerAmbiguous';
if isempty(rows)
    diagnostics = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    diagnostics = cell2table(rows, 'VariableNames', names);
end
end

function [epsT, epsAlternative, convention, alternativeConvention] = target_size(stage, cfg)
convention = string(cfg.equationOne.targetSizeConvention);
if strcmp(stage.targetKind, "final")
    width = rectangle_width(cfg.targetRectangle);
    height = cfg.targetRectangle.yMax - cfg.targetRectangle.yMin;
    diameterValue = min(width, height);
    radiusValue = 0.5 * diameterValue;
else
    diameterValue = 2 * cfg.proxyTargetRadius;
    radiusValue = cfg.proxyTargetRadius;
end
switch lower(convention)
    case "diameter"
        epsT = diameterValue;
        epsAlternative = radiusValue;
        alternativeConvention = "radius";
    case "radius"
        epsT = radiusValue;
        epsAlternative = diameterValue;
        alternativeConvention = "diameter";
    otherwise
        error('SchroerOtt:InvalidTargetSizeConvention', ...
            'cfg.equationOne.targetSizeConvention must be "diameter" or "radius".');
end
end

function w = rectangle_width(rect)
if rect.xMax >= rect.xMin
    w = rect.xMax - rect.xMin;
else
    w = rect.xMax + 1 - rect.xMin;
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
