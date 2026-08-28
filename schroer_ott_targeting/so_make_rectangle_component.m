function component = so_make_rectangle_component(rect, cfg, id)
%SO_MAKE_RECTANGLE_COMPONENT Closed rectangle boundary on S^1 x R.
if nargin < 3
    id = rect.id;
end
rect = apply_target_lift(rect, cfg);
component.type = 'rectangle';
component.id = char(id);
component.isClosed = true;
component.parameterRange = [0, 1];
width = rectangle_width(rect);
height = rect.yMax - rect.yMin;
component.characteristicSize = min(width, height);
component.rectangle = rect;
component.center = so_rectangle_center(rect);
component.chainID = "";
component.phasePointIndex = NaN;
component.gamma0 = @(s) rectangle_boundary(rect, s);
component.contains = @(z, tol) so_point_in_rectangle(z, rect, tol);
end

function rect = apply_target_lift(rect, cfg)
if isfield(rect, 'id') && strcmp(rect.id, 'final-target')
    rect.yMin = rect.yMin + cfg.transport.targetLiftShift;
    rect.yMax = rect.yMax + cfg.transport.targetLiftShift;
end
end

function w = rectangle_width(rect)
if rect.xMax >= rect.xMin
    w = rect.xMax - rect.xMin;
else
    w = rect.xMax + 1 - rect.xMin;
end
end

function z = rectangle_boundary(rect, s)
s = mod(reshape(s, 1, []), 1);
w = rectangle_width(rect);
h = rect.yMax - rect.yMin;
perim = 2 * (w + h);
ell = s * perim;
z = zeros(2, numel(s));
for i = 1:numel(s)
    a = ell(i);
    if a <= w
        x = rect.xMin + a;
        y = rect.yMin;
    elseif a <= w + h
        x = rect.xMin + w;
        y = rect.yMin + (a - w);
    elseif a <= 2 * w + h
        x = rect.xMin + w - (a - w - h);
        y = rect.yMax;
    else
        x = rect.xMin;
        y = rect.yMax - (a - 2 * w - h);
    end
    z(:, i) = [x; y];
end
end

