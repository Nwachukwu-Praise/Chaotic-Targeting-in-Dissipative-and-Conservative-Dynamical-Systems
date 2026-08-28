function result = shinbrot_parameter_targeting( ...
    initialState, target, params, control)
%SHINBROT_PARAMETER_TARGETING Target using one bounded parameter change.
%
% A search chooses a constant pControl in [-DeltaP, DeltaP]. That value is
% applied in the Lorenz differential equation from the source crossing
% until the target interval is reached. No state component is kicked.

search = search_parameter_to_target( ...
    initialState, target, params, control);

result.hit = false;
result.numCrossings = Inf;
result.physicalTime = Inf;
result.selectedParameter = search.selectedP;
result.controlAfterTarget = 0;
result.states = [];
result.xSequence = initialState(1);
result.eventTimes = [];
result.t = [];
result.X = [];
result.initialState = initialState(:);
result.search = search;

if ~search.found
    warning('No targeting parameter was found in the allowed interval.');
    result.search.controlledReplayCrossings = Inf;
    result.search.replayHit = false;
    result.search.replayFinalError = Inf;
    return;
end

states = NaN(3, search.horizon);
eventTimes = NaN(1, search.horizon);
tFull = [];
XFull = [];
currentState = initialState(:);
currentTime = 0;
numRecorded = 0;

for k = 1:search.horizon
    segment = next_valid_section_crossing( ...
        currentState, params, search.selectedP);

    if ~segment.success
        warning('Controlled replay failed before crossing %d.', k);
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
        result.hit = true;
        result.numCrossings = k;
        result.physicalTime = currentTime;
        break;
    end
end

states = states(:, 1:numRecorded);
eventTimes = eventTimes(1:numRecorded);
result.states = states;
result.xSequence = [initialState(1), states(1, :)];
result.eventTimes = eventTimes;
result.t = tFull;
result.X = XFull;
result.search.controlledReplayCrossings = result.numCrossings;
result.search.replayHit = result.hit;
if isempty(states)
    result.search.replayFinalError = Inf;
else
    result.search.replayFinalError = abs(states(1, end) - target.x);
end
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
