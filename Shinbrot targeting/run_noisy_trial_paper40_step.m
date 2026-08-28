function trial = run_noisy_trial_paper40_step( ...
    trial, sourceState, target, params, noiseControl, noise, ...
    pCurrent, verifiedIdentifier, trialTimer)
%RUN_NOISY_TRIAL_PAPER40_STEP Paper-mode noisy targeting trial.
%
% Retargeting is attempted after every 40 completed RK4 steps, using the
% current full state even when it is off the Poincare section.  Noise is
% applied only to the propagated plant, never inside deterministic forecast
% searches.

dt = noise.dt;
sigmaNoise = noise.sigmaNoise;
require_paper_retarget_interval(noise);
retargetStepInterval = noise.retargetStepInterval;
maxRk4Steps = require_field(noise, 'maxRk4Steps');
maxAcceptedCrossings = require_field(noise, 'maxControlledCrossings');
maxRejectedCrossings = require_field(noise, 'maxRejectedCrossings');
maxStateNorm = get_field_with_default(noise, 'maxStateNorm', Inf);

trial.retargetStepInterval = retargetStepInterval;
currentState = sourceState(:);
currentTime = 0;
tFull = 0;
XFull = currentState.';
states = NaN(3, maxAcceptedCrossings);
eventTimes = NaN(1, maxAcceptedCrossings);

for stepIndex = 1:maxRk4Steps
    hitEventState = [];
    hitEventTime = NaN;
    stepStartState = currentState + sigmaNoise * randn(3, 1);
    stepStartValue = stepStartState(3) - params.zSection;
    stepEndState = rk4_lorenz_step(stepStartState, dt, params, pCurrent);
    stepEndValue = stepEndState(3) - params.zSection;

    if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
        trial.integrationFailed = true;
        trial.failureReason = 'nonfinite or excessive state norm';
        break;
    end

    if crossed_in_direction(stepStartValue, stepEndValue, ...
            params.crossingDirection) && stepEndValue ~= stepStartValue
        alpha = -stepStartValue / (stepEndValue - stepStartValue);
        alpha = max(0, min(1, alpha));
        eventState = stepStartState + alpha * (stepEndState - stepStartState);
        eventState(3) = params.zSection;
        eventTime = currentTime + alpha * dt;

        if eventState(1) > params.xSectionMin
            trial.numAcceptedCrossings = trial.numAcceptedCrossings + 1;
            acceptedIndex = trial.numAcceptedCrossings;
            if acceptedIndex <= maxAcceptedCrossings
                states(:, acceptedIndex) = eventState;
                eventTimes(acceptedIndex) = eventTime;
            end
            trial.acceptedCrossingStates(:, end + 1) = eventState;
            trial.acceptedCrossingTimes(end + 1) = eventTime;
            trial.acceptedCrossingStepIndices(end + 1) = stepIndex;
            trial.acceptedCrossingX(end + 1) = eventState(1);
            trial.finalDistance = abs(eventState(1) - target.x);
            trial.finalState = eventState;

            if trial.finalDistance <= target.tolerance
                trial.hit = true;
                trial.numCrossingsToTarget = acceptedIndex;
                trial.numCrossings = acceptedIndex;
                trial.physicalTime = eventTime;
                hitEventState = eventState;
                hitEventTime = eventTime;
            end
        else
            trial.numRejectedCrossings = trial.numRejectedCrossings + 1;
            trial.rejectedCrossingStates(:, end + 1) = eventState;
            trial.rejectedCrossingTimes(end + 1) = eventTime;
            trial.rejectedCrossingStepIndices(end + 1) = stepIndex;
            if trial.numRejectedCrossings > maxRejectedCrossings
                trial.integrationFailed = true;
                trial.rejectedCrossingLimitExhausted = true;
                trial.failureReason = 'rejected-crossing limit exhausted';
            end
        end
    end

    currentTime = currentTime + dt;
    currentState = stepEndState;
    trial.numRk4Steps = stepIndex;
    tFull(end + 1, 1) = currentTime; %#ok<AGROW>
    XFull(end + 1, :) = currentState.'; %#ok<AGROW>

    if trial.hit || trial.integrationFailed
        if trial.hit && ~isempty(hitEventState)
            tFull(end) = hitEventTime;
            XFull(end, :) = hitEventState.';
        end
        break;
    end
    if trial.numAcceptedCrossings >= maxAcceptedCrossings
        trial.acceptedCrossingLimitExhausted = true;
        trial.failureReason = 'accepted-crossing limit exhausted';
        break;
    end

    if mod(stepIndex, retargetStepInterval) == 0
        newSearch = search_parameter_to_target( ...
            currentState, target, params, noiseControl, verifiedIdentifier);
        assert_verified_bisection_search( ...
            newSearch, verifiedIdentifier, ...
            'paper 40-step retargeting search', target);
        accepted = newSearch.found && isfinite(newSearch.selectedP);
        trial = record_noisy_search_event( ...
            trial, 'retarget', trial.numAcceptedCrossings, stepIndex, ...
            currentTime, newSearch, accepted);
        if accepted
            pCurrent = newSearch.selectedP;
            trial.selectedPValues(end + 1) = pCurrent;
        end
    end
end

if ~trial.hit && ~trial.integrationFailed && ...
        trial.numRk4Steps >= maxRk4Steps
    trial.stepBudgetExhausted = true;
    trial.failureReason = 'global RK4-step budget exhausted';
end

numStored = min(trial.numAcceptedCrossings, maxAcceptedCrossings);
trial.states = states(:, 1:numStored);
trial.eventTimes = eventTimes(1:numStored);
trial.xSequence = [sourceState(1), trial.states(1, :)];
trial.t = tFull;
trial.X = XFull;
trial.selectedPFinal = pCurrent;
trial.computationalTime = toc(trialTimer);
end

function require_paper_retarget_interval(noise)
if ~isfield(noise, 'retargetStepInterval') || ...
        noise.retargetStepInterval ~= 40
    error(['paper40StepRetargetMode requires ', ...
        'noise.retargetStepInterval == 40.']);
end
end

function value = require_field(s, fieldName)
if ~isfield(s, fieldName) || isempty(s.(fieldName)) || ...
        ~isfinite(s.(fieldName))
    error('Noise configuration must define finite %s.', fieldName);
end
value = s.(fieldName);
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
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
