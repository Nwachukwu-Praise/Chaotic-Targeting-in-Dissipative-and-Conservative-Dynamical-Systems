function component = so_make_control_component(zSource, cfg)
%SO_MAKE_CONTROL_COMPONENT Vertical admissible y-control segment.
delta = cfg.controlAmplitude;
component.type = 'control';
component.id = 'control-segment';
component.isClosed = false;
component.parameterRange = [-delta, delta];
component.characteristicSize = 2 * delta;
component.center = zSource(:);
component.phasePointIndex = NaN;
component.chainID = "";
component.gamma0 = @(u) [zSource(1) + zeros(1, numel(u)); zSource(2) + reshape(u, 1, [])];
component.contains = @(z, tol) abs(so_wrap_diff_x(z(1, :) - zSource(1))) <= tol & ...
    z(2, :) >= zSource(2) - delta - tol & z(2, :) <= zSource(2) + delta + tol;
end

