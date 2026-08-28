function segment = next_valid_section_crossing_noisy_rk4( ...
    X0, params, noise, pControl)
%NEXT_VALID_SECTION_CROSSING_NOISY_RK4 Find next noisy section crossing.
%
% Before each RK4 step, independent Gaussian coordinate noise is added:
%   stepStartState = currentState + sigmaNoise * randn(3, 1).
% The RK4 step is then taken from stepStartState to stepEndState. Section
% crossings are interpolated only across that integrated RK4 segment. An
% instantaneous crossing caused solely by the noise jump is deliberately not
% counted.

if nargin < 4
    pControl = 0;
end

dt = noise.dt;
if isfield(noise, 'sigmaNoise')
    sigmaNoise = noise.sigmaNoise;
elseif isfield(noise, 'noiseCoordinateSigma')
    sigmaNoise = noise.noiseCoordinateSigma;
elseif isfield(noise, 'delta')
    sigmaNoise = noise.delta;
else
    error('Noise configuration must define sigmaNoise.');
end
maxSteps = ceil(noise.maxTimeToNextValidSection / dt);

if isfield(noise, 'maxRejectedCrossings')
    maxRejectedCrossings = noise.maxRejectedCrossings;
else
    maxRejectedCrossings = 100;
end

if isfield(noise, 'maxStateNorm')
    maxStateNorm = noise.maxStateNorm;
else
    maxStateNorm = Inf;
end

currentState = X0(:);
currentTime = 0;

segment.success = false;
segment.failed = false;
segment.failureReason = '';
segment.t = 0;
segment.X = currentState.';
segment.eventTime = NaN;
segment.eventState = NaN(3, 1);
segment.x = NaN;
segment.numSteps = 0;
segment.numRejectedCrossings = 0;
segment.pControl = pControl;
segment.sigmaNoise = sigmaNoise;
segment.lastStepStartState = NaN(3, 1);
segment.lastStepEndState = NaN(3, 1);
segment.crossingUsesIntegratedSegmentEndpoints = true;
segment.instantaneousNoiseJumpCrossingsExcluded = true;

% Step clear of the section before testing for crossings.
%
% An accepted crossing is returned with its third coordinate set exactly to
% params.zSection, and the trial loop passes that state straight back in as
% the start of the next search.  Without a guard the first noise draw moves
% z off the section by order sigmaNoise and the first integrated step can
% carry it back across, so the crossing that was just returned is detected
% again a fraction of one step later.  The deterministic finder
% next_valid_section_crossing already avoids this by taking a short
% event-free step whenever the incoming state lies on the section; this is
% the fixed-step equivalent.
%
% Symptom without the guard: successive accepted crossings separated by far
% less than dt, against a mean return time of order one.
if isfield(params, 'sectionTolerance')
    sectionTolerance = params.sectionTolerance;
else
    sectionTolerance = 1e-9;
end
if isfield(noise, 'sectionClearanceSteps')
    clearanceSteps = noise.sectionClearanceSteps;
else
    % One step is enough to leave the section deterministically, but the
    % noise draw is applied before the step, so two steps are used to keep
    % the clearance large compared with a single sigmaNoise displacement.
    clearanceSteps = 2;
end
segment.sectionClearanceStepsUsed = 0;

if abs(currentState(3) - params.zSection) <= max(sectionTolerance, 0)
    for clearanceIndex = 1:clearanceSteps
        stepStartState = currentState + sigmaNoise * randn(3, 1);
        stepEndState = rk4_lorenz_step(stepStartState, dt, params, pControl);
        currentTime = currentTime + dt;
        segment.numSteps = segment.numSteps + 1;
        segment.sectionClearanceStepsUsed = clearanceIndex;
        segment.t(end + 1, 1) = currentTime;
        segment.X(end + 1, :) = stepEndState.';
        if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
            segment.failed = true;
            segment.failureReason = ...
                'state became nonfinite or too large while clearing the section';
            return;
        end
        currentState = stepEndState;
    end
end

for stepIndex = 1:maxSteps
    stepStartState = currentState + sigmaNoise * randn(3, 1);
    stepStartValue = stepStartState(3) - params.zSection;

    stepEndState = rk4_lorenz_step(stepStartState, dt, params, pControl);
    stepEndValue = stepEndState(3) - params.zSection;

    currentTime = currentTime + dt;
    segment.numSteps = segment.sectionClearanceStepsUsed + stepIndex;
    segment.lastStepStartState = stepStartState;
    segment.lastStepEndState = stepEndState;
    segment.t(end + 1, 1) = currentTime;
    segment.X(end + 1, :) = stepEndState.';

    if any(~isfinite(stepEndState)) || norm(stepEndState) > maxStateNorm
        segment.failed = true;
        segment.failureReason = 'state became nonfinite or too large';
        return;
    end

    if crossed_in_direction( ...
            stepStartValue, stepEndValue, params.crossingDirection)
        denominator = stepEndValue - stepStartValue;
        if denominator ~= 0
            alpha = -stepStartValue / denominator;
            alpha = max(0, min(1, alpha));
            eventState = stepStartState + ...
                alpha * (stepEndState - stepStartState);
            eventState(3) = params.zSection;
            eventTime = currentTime - dt + alpha * dt;

            if eventState(1) > params.xSectionMin
                segment.success = true;
                segment.eventTime = eventTime;
                segment.eventState = eventState;
                segment.x = eventState(1);
                segment.t(end) = eventTime;
                segment.X(end, :) = eventState.';
                return;
            end

            segment.numRejectedCrossings = ...
                segment.numRejectedCrossings + 1;
            if segment.numRejectedCrossings > maxRejectedCrossings
                segment.failed = true;
                segment.failureReason = ...
                    'too many rejected half-plane crossings';
                return;
            end
        end
    end

    currentState = stepEndState;
end

segment.failureReason = 'no valid section crossing before max time';
end

function tf = crossed_in_direction(valueOld, valueNew, direction)
%CROSSED_IN_DIRECTION True when a fixed-step pair crosses the section.

if direction > 0
    tf = (valueOld < 0) && (valueNew >= 0);
elseif direction < 0
    tf = (valueOld > 0) && (valueNew <= 0);
else
    tf = (valueOld == 0 && valueNew ~= 0) || ...
        (valueOld * valueNew < 0);
end
end
