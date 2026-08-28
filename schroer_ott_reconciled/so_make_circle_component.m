function component = so_make_circle_component(center, radius, id, chainID, phasePointIndex)
%SO_MAKE_CIRCLE_COMPONENT Circular target ball in cylinder metric.
component.type = 'circle';
component.id = char(id);
component.isClosed = true;
component.parameterRange = [0, 1];
component.characteristicSize = 2 * radius;
component.center = center(:);
component.radius = radius;
component.chainID = string(chainID);
component.phasePointIndex = phasePointIndex;
component.gamma0 = @(s) circle_boundary(center(:), radius, s);
component.contains = @(z, tol) so_cylinder_distance(z, center(:)) <= radius + tol;
end

function z = circle_boundary(center, radius, s)
s = mod(reshape(s, 1, []), 1);
theta = 2 * pi .* s;
z = [center(1) + radius .* cos(theta); center(2) + radius .* sin(theta)];
end

