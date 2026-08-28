function trial = run_noisy_trial_section_retarget( ...
    trial, sourceState, target, params, noiseControl, noise, ...
    pCurrent, verifiedIdentifier, trialTimer)
%RUN_NOISY_TRIAL_SECTION_RETARGET Noisy section-based retargeting extension.
%
% This is not the paper's stated 40-step retargeting rule.  It is retained
% as a project extension: after an accepted Poincare crossing misses the
% target, pControl is recomputed from that section state.

dt = noise.dt;
sigmaNoise = noise.sigmaNoise;
maxRk4Steps = require_field(noise, 'maxRk4Steps');
maxAcceptedCrossings = require_field(noise, 'maxControlledCrossings');
maxRejectedCrossings = require_field(noise, 'maxRejectedCrossings');
maxStateNorm = get_field_with_default(noise, 'maxStateNorm', Inf);
sectionClearanceSteps = get_field_with_default(noise, ...
    'sectionClearanceSteps', 2);
sectionTolerance = get_field_with_default(params, ...
    'sectionTolerance', 1e-9);

currentState = sourceState(:);
currentTime = 0;
tFull = 0;
XFull = currentState.';
states = NaN(3, maxAcceptedCrossings);
eventTimes = NaN(1, maxAcceptedCrossings);

while trial.numRk4Steps < maxRk4Steps
    if abs(currentState(3) - params.zSection) <= sectionTolerance
        [currentState, currentTime, tFull, XFull, trial] = ...
            advance_without_event_detection(currentState, currentTime, ...
            tFull, XFull, trial, params, noise, pCurrent, ...
            maxRk4Steps, sectionClearanceSteps, maxStateNorm);
        if trial.integrationFailed || trial.numRk4Steps >= maxRk4Steps
            break;
        end
    end

    stepStartState = currentState + sigmaNoise * randn(3, 1);
    stepStartValue = stepStartState(3) - params.zSection;
    stepEndState = rk4_lorenz_step(stepStartState, dt, params, pCurrent);
    stepEndValue = stepEndState(3) - params.zSection;

    if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
        trial.integrationFailed = true;
        trial.failureReason = 'nonfinite or excessive state norm';
        break;
    end

    eventDetected = crossed_in_direction( ...
        stepStartValue, stepEndValue, params.crossingDirection) && ...
        stepEndValue ~= stepStartValue;

    trial.numRk4Steps = trial.numRk4Steps + 1;
    currentTime = currentTime + dt;
    currentState = stepEndState;
    tFull(end + 1, 1) = currentTime; %#ok<AGROW>
    XFull(end + 1, :) = currentState.'; %#ok<AGROW>

    if eventDetected
        alpha = -stepStartValue / (stepEndValue - stepStartValue);
        alpha = max(0, min(1, alpha));
        eventState = stepStartState + alpha * (stepEndState - stepStartState);
        eventState(3) = params.zSection;
        eventTime = currentTime - dt + alpha * dt;

        if eventState(1) > params.xSectionMin
            trial.numAcceptedCrossings = trial.numAcceptedCrossings + 1;
            acceptedIndex = trial.numAcceptedCrossings;
            if acceptedIndex <= maxAcceptedCrossings
                states(:, acceptedIndex) = eventState;
                eventTimes(acceptedIndex) = eventTime;
            end
            trial.acceptedCrossingStates(:, end + 1) = eventState;
            trial.acceptedCrossingTimes(end + 1) = eventTime;
            trial.acceptedCrossingStepIndices(end + 1) = trial.numRk4Steps;
            trial.acceptedCrossingX(end + 1) = eventState(1);
            trial.finalDistance = abs(eventState(1) - target.x);
            trial.finalState = eventState;

            tFull(end) = eventTime;
            XFull(end, :) = eventState.';
            currentTime = eventTime;
            currentState = eventState;

            if trial.finalDistance <= target.tolerance
                trial.hit = true;
                trial.numCrossingsToTarget = acceptedIndex;
                trial.numCrossings = acceptedIndex;
                trial.physicalTime = eventTime;
                break;
            end

            if acceptedIndex >= maxAcceptedCrossings
                trial.acceptedCrossingLimitExhausted = true;
                trial.failureReason = 'accepted-crossing limit exhausted';
                break;
            end

            newSearch = search_parameter_to_target( ...
                currentState, target, params, noiseControl, ...
                verifiedIdentifier);
            assert_verified_bisection_search(newSearch, verifiedIdentifier, ...
                'section retargeting search', target);
            accepted = newSearch.found && isfinite(newSearch.selectedP);
            trial = record_noisy_search_event( ...
                trial, 'retarget', acceptedIndex, trial.numRk4Steps, ...
                currentTime, newSearch, accepted);
            if accepted
                pCurrent = newSearch.selectedP;
                trial.selectedPValues(end + 1) = pCurrent;
            end
        else
            trial.numRejectedCrossings = trial.numRejectedCrossings + 1;
            trial.rejectedCrossingStates(:, end + 1) = eventState;
            trial.rejectedCrossingTimes(end + 1) = eventTime;
            trial.rejectedCrossingStepIndices(end + 1) = trial.numRk4Steps;
            if trial.numRejectedCrossings > maxRejectedCrossings
                trial.integrationFailed = true;
                trial.rejectedCrossingLimitExhausted = true;
                trial.failureReason = 'rejected-crossing limit exhausted';
                break;
            end
        end
    end
end

if ~trial.hit && ~trial.integrationFailed && ...
        trial.numRk4Steps >= maxRk4Steps && ...
        isempty(trial.failureReason)
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

function [currentState, currentTime, tFull, XFull, trial] = ...
    advance_without_event_detection(currentState, currentTime, tFull, ...
    XFull, trial, params, noise, pCurrent, maxRk4Steps, ...
    clearanceSteps, maxStateNorm)
for k = 1:clearanceSteps
    if trial.numRk4Steps >= maxRk4Steps
        return;
    end
    stepStartState = currentState + noise.sigmaNoise * randn(3, 1);
    stepEndState = rk4_lorenz_step(stepStartState, noise.dt, params, pCurrent);
    if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
        trial.integrationFailed = true;
        trial.failureReason = ...
            'nonfinite or excessive state norm while clearing section';
        return;
    end
    trial.numRk4Steps = trial.numRk4Steps + 1;
    currentTime = currentTime + noise.dt;
    currentState = stepEndState;
    tFull(end + 1, 1) = currentTime; %#ok<AGROW>
    XFull(end + 1, :) = currentState.'; %#ok<AGROW>
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
