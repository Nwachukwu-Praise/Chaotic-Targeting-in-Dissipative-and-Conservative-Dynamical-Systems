function du = lorenz_ode(u, p)
%LORENZ_ODE  Right hand side of the Lorenz system (Table 2 of the paper).
%
%   du = LORENZ_ODE(U,P) with U = [x;y;z] and P = [sigma;rho;beta]:
%
%       dx/dt = sigma*(y - x)
%       dy/dt = rho*x - y - x*z
%       dz/dt = x*y - beta*z
%
%   Nominal parameters: sigma = 10, rho = 28, beta = 8/3.

du = [ p(1)*(u(2) - u(1)); ...
       p(2)*u(1) - u(2) - u(1)*u(3); ...
       u(1)*u(2) - p(3)*u(3) ];
end
