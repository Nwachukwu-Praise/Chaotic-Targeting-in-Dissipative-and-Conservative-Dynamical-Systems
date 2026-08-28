function replay = so_replay_with_noise(result, sigma, noiseSequence)
%SO_REPLAY_WITH_NOISE Re-run a fixed control schedule under additive noise.
%
%   replay = SO_REPLAY_WITH_NOISE(result, sigma)
%   replay = SO_REPLAY_WITH_NOISE(result, sigma, eta)
%
% The control values and the iteration counts are taken from the
% deterministic solution and NOT recomputed.  That is deliberate: it is the
% open-loop test, and it is the strict one.  Schroer and Ott's claim is that
% their trajectories survive small noise without any change to the scheme
% ("Without any change to our targeting scheme, all trajectories reached the
% target"), and for larger noise they repeat the forward-backward procedure
% every few steps.  Replaying the schedule unchanged measures the first
% claim; it will fail before a closed-loop scheme would, which is the point.
%
% Noise enters the momentum update, the same place the control does:
%     y_{n+1} = y_n - (k/2pi) sin(2 pi x_n) + eta_n
%     x_{n+1} = x_n + y_{n+1}
%
% eta may be supplied so that several sigma levels can share one draw.

cfg = result.configuration;
segments = result.executionSegments;
totalIterations = sum([segments.iterations]);
if nargin < 3 || isempty(noiseSequence)
    noiseSequence = randn(1, totalIterations);
end
if numel(noiseSequence) < totalIterations
    error('SchroerOtt:ShortNoiseSequence', ...
        'Need %d noise samples, received %d.', totalIterations, numel(noiseSequence));
end

z = result.sourceState;
trajectory = z;
kicks = zeros(2, 0);
noiseIndex = 0;
for i = 1:numel(segments)
    z = z + [0; segments(i).control];
    kicks(:, end + 1) = z; %#ok<AGROW>
    trajectory = [trajectory, [NaN; NaN], z]; %#ok<AGROW>
    for n = 1:segments(i).iterations
        noiseIndex = noiseIndex + 1;
        eta = sigma * noiseSequence(noiseIndex);
        yNext = z(2) - cfg.k / (2 * pi) * sin(2 * pi * z(1)) + eta;
        z = [z(1) + yNext; yNext];
        trajectory(:, end + 1) = z; %#ok<AGROW>
    end
end

targetComponent = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
contained = targetComponent.contains(z, cfg.containmentTolerance);
deterministicEnd = result.finalState;

replay.sigma = sigma;
replay.totalIterations = totalIterations;
replay.noiseSequence = noiseSequence(1:totalIterations);
replay.trajectory = trajectory;
replay.controlKickStates = kicks;
replay.finalState = z;
replay.targetContained = contained;
replay.endpointDisplacement = so_cylinder_distance(z, deterministicEnd);
end
