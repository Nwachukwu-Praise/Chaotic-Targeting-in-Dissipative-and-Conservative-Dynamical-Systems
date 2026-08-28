function data = generate_return_map(X0, numCrossings, params, pControl)
%GENERATE_RETURN_MAP Generate accepted crossings and the scalar return map.
%
% Full crossing states are retained, while the approximately one-
% dimensional map is represented by consecutive x coordinates:
%       X_{n+1} = F(pControl, X_n).

if nargin < 4
    pControl = 0;
end

states = NaN(3, numCrossings);
eventTimes = NaN(1, numCrossings);
tFull = [];
XFull = [];

currentState = X0(:);
currentTime = 0;
numGenerated = 0;

for k = 1:numCrossings
    segment = next_valid_section_crossing( ...
        currentState, params, pControl);

    if ~segment.success
        warning('Return map stopped before accepted crossing %d.', k);
        break;
    end

    [tFull, XFull] = append_segment( ...
        tFull, XFull, currentTime, segment);

    currentTime = currentTime + segment.eventTime;
    currentState = segment.eventState;
    numGenerated = k;

    states(:, k) = currentState;
    eventTimes(k) = currentTime;
end

states = states(:, 1:numGenerated);
eventTimes = eventTimes(1:numGenerated);

data.states = states;
data.x = states(1, :);
if numGenerated >= 2
    data.xN = data.x(1:end-1);
    data.xNext = data.x(2:end);
else
    data.xN = [];
    data.xNext = [];
end
data.eventTimes = eventTimes;
data.t = tFull;
data.X = XFull;
data.pControl = pControl;
data.initialState = X0(:);
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
