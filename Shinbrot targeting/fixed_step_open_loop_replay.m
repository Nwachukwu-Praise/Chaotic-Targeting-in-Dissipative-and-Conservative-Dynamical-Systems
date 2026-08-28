function comparison = fixed_step_open_loop_replay( ...
    sourceState, target, params, noise, pControl, horizon)
%FIXED_STEP_OPEN_LOOP_REPLAY Compare adaptive and fixed-step zero-noise maps.
%
% The same constant pControl is used in both calculations.  No retargeting
% is performed.

adaptive = adaptive_replay(sourceState, target, params, pControl, horizon);
fixed = fixed_rk4_replay(sourceState, target, params, noise, ...
    pControl, horizon);

comparison.selectedP = pControl;
comparison.horizon = horizon;
comparison.adaptiveFinalX = adaptive.finalX;
comparison.fixedStepFinalX = fixed.finalX;
comparison.adaptiveTargetError = adaptive.targetError;
comparison.fixedStepTargetError = fixed.targetError;
comparison.absoluteSolverDiscrepancy = ...
    abs(fixed.finalX - adaptive.finalX);
comparison.discrepancyFractionOfTolerance = ...
    comparison.absoluteSolverDiscrepancy / target.tolerance;
comparison.rk4StepCount = fixed.numRk4Steps;
comparison.acceptedCrossingCount = fixed.numAcceptedCrossings;
comparison.rejectedCrossingCount = fixed.numRejectedCrossings;
comparison.adaptiveHit = adaptive.hit;
comparison.fixedStepHit = fixed.hit;
comparison.fixedStepFailureReason = fixed.failureReason;
comparison.fixedStepFinalState = fixed.finalState;
comparison.adaptiveFinalState = adaptive.finalState;
tableStruct = rmfield(comparison, ...
    {'fixedStepFinalState', 'adaptiveFinalState'});
comparison.table = struct2table(tableStruct, 'AsArray', true);
end

function replay = adaptive_replay(sourceState, target, params, pControl, horizon)
currentState = sourceState(:);
replay.finalState = NaN(3, 1);
replay.finalX = NaN;
replay.targetError = NaN;
replay.hit = false;
for k = 1:horizon
    segment = next_valid_section_crossing(currentState, params, pControl);
    if ~segment.success
        return;
    end
    currentState = segment.eventState;
end
replay.finalState = currentState;
replay.finalX = currentState(1);
replay.targetError = abs(currentState(1) - target.x);
replay.hit = replay.targetError <= target.tolerance;
end

function replay = fixed_rk4_replay( ...
    sourceState, target, params, noise, pControl, horizon)
dt = noise.dt;
maxRk4Steps = get_field_with_default(noise, 'maxRk4Steps', Inf);
maxRejectedCrossings = get_field_with_default(noise, ...
    'maxRejectedCrossings', Inf);
maxStateNorm = get_field_with_default(noise, 'maxStateNorm', Inf);

currentState = sourceState(:);
currentTime = 0; %#ok<NASGU>
numAccepted = 0;
numRejected = 0;
replay.numRk4Steps = 0;
replay.finalState = currentState;
replay.finalX = NaN;
replay.targetError = NaN;
replay.hit = false;
replay.failureReason = '';

while replay.numRk4Steps < maxRk4Steps && numAccepted < horizon
    stepStartState = currentState;
    stepStartValue = stepStartState(3) - params.zSection;
    stepEndState = rk4_lorenz_step(stepStartState, dt, params, pControl);
    stepEndValue = stepEndState(3) - params.zSection;
    replay.numRk4Steps = replay.numRk4Steps + 1;

    if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
        replay.failureReason = 'nonfinite or excessive state norm';
        break;
    end

    if crossed_in_direction(stepStartValue, stepEndValue, ...
            params.crossingDirection) && stepEndValue ~= stepStartValue
        alpha = -stepStartValue / (stepEndValue - stepStartValue);
        alpha = max(0, min(1, alpha));
        eventState = stepStartState + alpha * (stepEndState - stepStartState);
        eventState(3) = params.zSection;
        if eventState(1) > params.xSectionMin
            numAccepted = numAccepted + 1;
            currentState = eventState;
            replay.finalState = eventState;
            continue;
        end
        numRejected = numRejected + 1;
        if numRejected > maxRejectedCrossings
            replay.failureReason = 'rejected-crossing limit exhausted';
            break;
        end
    end

    currentState = stepEndState;
end

replay.numAcceptedCrossings = numAccepted;
replay.numRejectedCrossings = numRejected;
if numAccepted >= horizon
    replay.finalX = replay.finalState(1);
    replay.targetError = abs(replay.finalX - target.x);
    replay.hit = replay.targetError <= target.tolerance;
elseif isempty(replay.failureReason)
    replay.failureReason = 'global RK4-step budget exhausted';
end
end

function tf = crossed_in_direction(valueOld, valueNew, direction)
if direction > 0
    tf = (valueOld < 0) && (valueNew >= 0);
elseif direction < 0
    tf = (valueOld > 0) && (valueNew <= 0);
else
    tf = (valueOld == 0 && valueNew ~= 0) || ...
        (valueOld * valueNew < 0);
end
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
