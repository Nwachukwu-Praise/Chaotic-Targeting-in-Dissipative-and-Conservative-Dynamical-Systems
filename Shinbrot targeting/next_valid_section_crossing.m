function segment = next_valid_section_crossing(X0, params, pControl)
%NEXT_VALID_SECTION_CROSSING Advance to the next accepted Poincare crossing.
%
% A valid map iterate crosses z = params.zSection in the selected direction
% and satisfies x > params.xSectionMin. Plane crossings with smaller x are
% rejected, after which integration continues to the next plane crossing.
%
% If X0 lies on the section, a tiny event-free step is taken first. This
% prevents ode45 from immediately rediscovering the starting crossing.

if nargin < 3
    pControl = 0;
end

segment.success = false;
segment.t = [];
segment.X = [];
segment.eventTime = NaN;
segment.eventState = NaN(3, 1);
segment.x = NaN;
segment.numRejectedCrossings = 0;

currentState = X0(:);
currentTime = 0;

eventOptions = odeset(params.odeOptions, ...
    'Events', @(t, X) poincare_event(t, X, params));
noEventOptions = odeset(params.odeOptions, ...
    'Events', [], 'MaxStep', params.maxStepOffSection);

for crossingAttempt = 1:(params.maxRejectedCrossings + 1)
    if abs(currentState(3) - params.zSection) <= params.sectionTolerance
        stepEnd = min( ...
            currentTime + params.eventDisableTime, ...
            params.maxTimeToNextValidSection);

        if stepEnd <= currentTime
            return;
        end

        [tStep, XStep] = ode45( ...
            @(t, X) lorenz_rhs(t, X, params, pControl), ...
            [currentTime, stepEnd], currentState, noEventOptions);

        [segment.t, segment.X] = append_samples( ...
            segment.t, segment.X, tStep, XStep);
        currentTime = tStep(end);
        currentState = XStep(end, :).';
    end

    if currentTime >= params.maxTimeToNextValidSection
        return;
    end

    [tEvent, XEvent, te, Xe] = ode45( ...
        @(t, X) lorenz_rhs(t, X, params, pControl), ...
        [currentTime, params.maxTimeToNextValidSection], ...
        currentState, eventOptions);

    [segment.t, segment.X] = append_samples( ...
        segment.t, segment.X, tEvent, XEvent);

    if isempty(te)
        return;
    end

    currentTime = te(end);
    currentState = Xe(end, :).';

    if currentState(1) > params.xSectionMin
        segment.success = true;
        segment.eventTime = currentTime;
        segment.eventState = currentState;
        segment.x = currentState(1);
        return;
    end

    segment.numRejectedCrossings = segment.numRejectedCrossings + 1;
end
end

function [tFull, XFull] = append_samples(tFull, XFull, tNew, XNew)
%APPEND_SAMPLES Join ode45 output without duplicating boundary samples.

if isempty(tFull)
    tFull = tNew;
    XFull = XNew;
else
    tFull = [tFull; tNew(2:end)];
    XFull = [XFull; XNew(2:end, :)];
end
end
