function z = so_curve_point(component, parameter, direction, iterateCount, cfg)
%SO_CURVE_POINT Evaluate F^n(gamma0(s)) or F^{-n}(gamma0(s)).
z0 = component.gamma0(parameter);
z = so_iterate(z0, iterateCount, cfg, direction, false);
end

