function segment = integrate_noisy_lorenz_fixed_step( ...
    X0, numSteps, dt, params, pControl, sigmaNoise)
%INTEGRATE_NOISY_LORENZ_FIXED_STEP Fixed-step noisy Lorenz integration.
%
% Before every RK4 step, independent Gaussian noise with standard deviation
% sigmaNoise is added directly to X, Y and Z:
%   stepStartState = currentState + sigmaNoise * randn(3, 1).
% The RK4 step is integrated from stepStartState. This helper deliberately
% avoids ode45 because adaptive internal steps do not give a clean meaning
% to per-step coordinate noise.

if nargin < 5
    pControl = 0;
end
if nargin < 6
    sigmaNoise = 0;
end

X = X0(:);
t = (0:numSteps).' * dt;
states = NaN(numSteps + 1, 3);
states(1, :) = X.';
failed = false;
failureReason = '';

for k = 1:numSteps
    stepStartState = X + sigmaNoise * randn(3, 1);
    X = rk4_lorenz_step(stepStartState, dt, params, pControl);

    if any(~isfinite(X))
        failed = true;
        failureReason = 'nonfinite state';
        states = states(1:k, :);
        t = t(1:k);
        break;
    end

    states(k + 1, :) = X.';
end

segment.t = t;
segment.X = states;
segment.finalState = states(end, :).';
segment.numSteps = size(states, 1) - 1;
segment.failed = failed;
segment.failureReason = failureReason;
segment.pControl = pControl;
segment.sigmaNoise = sigmaNoise;
end
