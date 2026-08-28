function U = flow_rk4_step_vec(odefun, U, P, h)
%FLOW_RK4_STEP_VEC  Vectorised RK4 step with per-column step size.
%
%   U = FLOW_RK4_STEP_VEC(ODEFUN,U,P,H) advances every column of U (nu x M)
%   by its own step H (scalar or 1 x M).  ODEFUN must be vectorised, i.e.
%   accept (nu x M, r x M) and return nu x M.

if isscalar(h)
    hh = h;
else
    hh = repmat(h(:).', size(U,1), 1);
end
k1 = odefun(U, P);
k2 = odefun(U + (hh/2).*k1, P);
k3 = odefun(U + (hh/2).*k2, P);
k4 = odefun(U + hh.*k3,     P);
U  = U + (hh/6).*(k1 + 2*k2 + 2*k3 + k4);
end
