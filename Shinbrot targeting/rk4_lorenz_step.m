function XNext = rk4_lorenz_step(X, dt, params, pControl)
%RK4_LORENZ_STEP One fourth-order Runge-Kutta step for the Lorenz system.
%
% The control parameter enters only through the second equation:
%   dy/dt = x*(rho - z) - y + pControl.

if nargin < 4
    pControl = 0;
end

X = X(:);
k1 = lorenz_rhs(0, X, params, pControl);
k2 = lorenz_rhs(0, X + 0.5 * dt * k1, params, pControl);
k3 = lorenz_rhs(0, X + 0.5 * dt * k2, params, pControl);
k4 = lorenz_rhs(0, X + dt * k3, params, pControl);

XNext = X + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end
