function u = flow_rk4_step(odefun, u, p, dt)
%FLOW_RK4_STEP  One classical Runge-Kutta 4 step of du/dt = odefun(u,p).
%
%   u = FLOW_RK4_STEP(ODEFUN, U, P, DT)
%
%   ODEFUN must accept (u,p) with u an nu x 1 column vector and return an
%   nu x 1 column vector.  Autonomous form is used throughout; time-dependent
%   forcing is handled by carrying the phase as an extra state (dz/dt = 1).

k1 = odefun(u,            p);
k2 = odefun(u + dt/2*k1,  p);
k3 = odefun(u + dt/2*k2,  p);
k4 = odefun(u + dt  *k3,  p);
u  = u + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
