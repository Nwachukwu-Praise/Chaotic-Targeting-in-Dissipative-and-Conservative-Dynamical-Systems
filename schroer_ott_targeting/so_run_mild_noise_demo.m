function noiseDemo = so_run_mild_noise_demo(deterministicResult, saveOutputs)
%SO_RUN_MILD_NOISE_DEMO Replay the deterministic schedule with one noise draw.
if nargin < 1 || isempty(deterministicResult)
    cfg = schroer_ott_default_config();
    cfg.saveFigures = false;
    catalogue = so_enumerate_periodic_orbits(cfg);
    route = so_construct_first_light_route(catalogue, cfg);
    deterministicResult = so_multistage_targeting(cfg, catalogue, route);
end
if nargin < 2
    saveOutputs = true;
end
cfg = deterministicResult.configuration;
delta = 2 * cfg.controlAmplitude;
sigma = 1e-2 * delta;
seed = 314159;
rng(seed);
totalIterations = sum([deterministicResult.executionSegments.iterations]);
noiseSequence = sigma .* randn(totalIterations, 1);

timer = tic;
[trajectory, finalState, replayError] = replay_with_noise(deterministicResult, cfg, noiseSequence);
targetComponent = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
targetContained = targetComponent.contains(finalState, cfg.containmentTolerance);
deterministicEndpoint = deterministicResult.finalState;
endpointDisplacement = norm([so_wrap_diff_x(finalState(1) - deterministicEndpoint(1)); ...
    finalState(2) - deterministicEndpoint(2)]);

noiseDemo.sigma = sigma;
noiseDemo.delta = delta;
noiseDemo.seed = seed;
noiseDemo.noiseApplication = "eta_n added to y_{n+1} before x_{n+1}=x_n+y_{n+1}";
noiseDemo.noiseSequence = noiseSequence;
noiseDemo.deterministicSuccess = deterministicResult.targetReached;
noiseDemo.noisySuccess = targetContained;
noiseDemo.totalIterations = totalIterations;
noiseDemo.maximumControl = max(abs(deterministicResult.executedControls.controlY));
noiseDemo.storedNoiseReplayError = replayError;
noiseDemo.finalState = finalState;
noiseDemo.finalContainment = targetContained;
noiseDemo.endpointDisplacement = endpointDisplacement;
noiseDemo.trajectory = trajectory;
noiseDemo.runtimeSeconds = toc(timer);
if saveOutputs
    if ~exist(cfg.outputDirectory, 'dir')
        mkdir(cfg.outputDirectory);
    end
    save(fullfile(cfg.outputDirectory, 'mild_noise_demo.mat'), 'noiseDemo', '-v7.3');
end
end

function [trajectory, finalState, replayError] = replay_with_noise(result, cfg, noiseSequence)
segments = result.executionSegments;
z = so_rectangle_center(cfg.sourceRectangle);
trajectory = z;
noiseIndex = 0;
for i = 1:numel(segments)
    z = z + [0; segments(i).control];
    trajectory = [trajectory, [NaN; NaN], z]; %#ok<AGROW>
    for n = 1:segments(i).iterations
        noiseIndex = noiseIndex + 1;
        z = noisy_step(z, cfg, noiseSequence(noiseIndex));
        trajectory(:, end + 1) = z; %#ok<AGROW>
    end
end
finalState = z;
trajectoryCheck = replay_again(result, cfg, noiseSequence);
d = trajectory - trajectoryCheck;
finiteColumns = all(isfinite(d), 1);
replayError = max(vecnorm(d(:, finiteColumns)));
end

function zNext = noisy_step(z, cfg, eta)
yNext = z(2) - cfg.k / (2 * pi) * sin(2 * pi * z(1)) + eta;
xNext = z(1) + yNext;
zNext = [xNext; yNext];
end

function trajectory = replay_again(result, cfg, noiseSequence)
segments = result.executionSegments;
z = so_rectangle_center(cfg.sourceRectangle);
trajectory = z;
noiseIndex = 0;
for i = 1:numel(segments)
    z = z + [0; segments(i).control];
    trajectory = [trajectory, [NaN; NaN], z]; %#ok<AGROW>
    for n = 1:segments(i).iterations
        noiseIndex = noiseIndex + 1;
        z = noisy_step(z, cfg, noiseSequence(noiseIndex));
        trajectory(:, end + 1) = z; %#ok<AGROW>
    end
end
end
