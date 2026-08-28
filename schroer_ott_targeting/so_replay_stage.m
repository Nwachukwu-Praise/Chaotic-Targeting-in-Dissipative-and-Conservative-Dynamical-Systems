function replay = so_replay_stage(zPre, controlY, plannedPath, selectedIterations, cfg, targetComponents)
%SO_REPLAY_STAGE Independently replay one accepted controlled stage.
%
% Replay deliberately does not call the connection solver and does not use
% so_iterate.  It applies the stored control once, advances the lifted
% standard map one step at a time, and compares every stored planned state
% with the replayed state.
if nargin < 6
    targetComponents = {};
end
if isempty(plannedPath)
    plannedPath = NaN(2, 0);
end

timer = tic;
zPre = zPre(:);
zPost = zPre + [0; controlY];
selectedIterations = max(0, selectedIterations);
replayPath = zeros(2, selectedIterations + 1);
replayPath(:, 1) = zPost;
z = zPost;
for n = 1:selectedIterations
    z = so_standard_map_lifted(z, cfg);
    replayPath(:, n + 1) = z;
end

storedCount = min(size(plannedPath, 2), selectedIterations + 1);
pathErrors = Inf(1, selectedIterations + 1);
if storedCount > 0
    for i = 1:storedCount
        pathErrors(i) = cylindrical_state_error(replayPath(:, i), plannedPath(:, i));
    end
end
if storedCount < selectedIterations + 1
    pathErrors(storedCount + 1:end) = Inf;
end

switchIndex = min(selectedIterations + 1, size(replayPath, 2));
switchStateReplayError = pathErrors(switchIndex);
finalStateReplayError = pathErrors(end);
targetContained = contains_any(targetComponents, replayPath(:, end), cfg);
controlWithinBound = abs(controlY) <= cfg.controlAmplitude + 10 * eps(cfg.controlAmplitude);

replay.preControlState = zPre;
replay.controlY = controlY;
replay.postControlState = zPost;
replay.selectedIterations = selectedIterations;
replay.replayPath = replayPath;
replay.storedPlannedPath = plannedPath(:, 1:min(size(plannedPath, 2), selectedIterations + 1));
replay.pathwiseReplayErrors = pathErrors;
replay.maxPathwiseReplayError = max(pathErrors);
replay.switchState = replayPath(:, switchIndex);
replay.finalState = replayPath(:, end);
replay.switchStateReplayError = switchStateReplayError;
replay.finalStateReplayError = finalStateReplayError;
replay.targetContained = targetContained;
replay.controlWithinBound = controlWithinBound;
replay.replayPassed = controlWithinBound && replay.maxPathwiseReplayError <= cfg.propagationConsistencyTolerance;
replay.runtimeSeconds = toc(timer);
end

function err = cylindrical_state_error(a, b)
err = norm([so_wrap_diff_x(a(1) - b(1)); a(2) - b(2)]);
end

function tf = contains_any(components, z, cfg)
tf = false;
for i = 1:numel(components)
    if components{i}.contains(z, cfg.containmentTolerance)
        tf = true;
        return;
    end
end
end
