function dU = lorenz_ode_vec(U, P)
%LORENZ_ODE_VEC  Vectorised Lorenz right hand side.
%
%   dU = LORENZ_ODE_VEC(U,P) with U (3 x M) = [x;y;z] and P (3 x M) =
%   [sigma;rho;beta] evaluates M copies of the Lorenz system at once.  All
%   copies are integrated in lock-step by FLOW_POINCARE_STEP_VEC, which is
%   what makes the Poincare map fast enough for the ECR experiments.

dU = [ P(1,:).*(U(2,:) - U(1,:)); ...
       P(2,:).*U(1,:) - U(2,:) - U(1,:).*U(3,:); ...
       U(1,:).*U(2,:) - P(3,:).*U(3,:) ];
end
