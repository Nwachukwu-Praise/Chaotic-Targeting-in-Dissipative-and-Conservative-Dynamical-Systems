function result = uncontrolled_hitting_time_1d( ...
    initialState, target, params, maxCrossings)
%UNCONTROLLED_HITTING_TIME_1D Measure natural recurrence to the target.
%
% Starting from the same source used for control, this baseline uses
% pControl = 0 and counts accepted Poincare crossings until
% abs(X_n - X_t) <= target.tolerance.

states = NaN(3, maxCrossings);
eventTimes = NaN(1, maxCrossings);
tFull = [];
XFull = [];

currentState = initialState(:);
currentTime = 0;
numRecorded = 0;
hitIndex = Inf;

for k = 1:maxCrossings
    segment = next_valid_section_crossing(currentState, params, 0);

    if ~segment.success
        warning('Uncontrolled run failed before crossing %d.', k);
        break;
    end

    [tFull, XFull] = append_segment( ...
        tFull, XFull, currentTime, segment);

    currentTime = currentTime + segment.eventTime;
    currentState = segment.eventState;
    numRecorded = k;

    states(:, k) = currentState;
    eventTimes(k) = currentTime;

    if abs(currentState(1) - target.x) <= target.tolerance
        hitIndex = k;
        break;
    end
end

states = states(:, 1:numRecorded);
eventTimes = eventTimes(1:numRecorded);

result.hit = isfinite(hitIndex);
result.numCrossings = hitIndex;
if result.hit
    result.physicalTime = eventTimes(hitIndex);
else
    result.numCrossings = Inf;
    result.physicalTime = Inf;
end
result.states = states;
result.xSequence = [initialState(1), states(1, :)];
result.eventTimes = eventTimes;
result.t = tFull;
result.X = XFull;
result.initialState = initialState(:);
end

function [tFull, XFull] = append_segment( ...
    tFull, XFull, currentTime, segment)
%APPEND_SEGMENT Convert segment-local times to one continuous time array.

tSegment = currentTime + segment.t;

if isempty(tFull)
    tFull = tSegment;
    XFull = segment.X;
else
    tFull = [tFull; tSegment(2:end)];
    XFull = [XFull; segment.X(2:end, :)];
end
end
