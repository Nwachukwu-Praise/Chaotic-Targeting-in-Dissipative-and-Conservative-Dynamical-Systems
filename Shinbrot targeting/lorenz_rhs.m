function dXdt = lorenz_rhs(~, X, params, pControl)
%LORENZ_RHS Lorenz equations with the Shinbrot control parameter.
%
% Control acts only through the additive parameter pControl in dy/dt.
% The state itself is never kicked or reset.

if nargin < 4
    pControl = 0;
end

x = X(1);
y = X(2);
z = X(3);

dXdt = [
    params.sigma * (y - x);
    x * (params.rho - z) - y + pControl;
    x * y - params.beta * z
    ];
end
