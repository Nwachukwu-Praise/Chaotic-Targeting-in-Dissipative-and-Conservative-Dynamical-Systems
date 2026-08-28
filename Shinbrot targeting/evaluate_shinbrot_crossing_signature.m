function result = evaluate_shinbrot_crossing_signature( ...
    sourceState, pControl, horizon, params)
%EVALUATE_SHINBROT_CROSSING_SIGNATURE Evaluate X_n(p) and crossing order.
%
% This helper is dedicated to the discontinuity-aware Shinbrot benchmark.
% Unlike next_valid_section_crossing, it records every directed crossing of
% the section plane. A left crossing (x < xSectionMin) is stored as -1 and
% an accepted right crossing (x > xSectionMin) is stored as +1.

result.xFinal = NaN;
result.success = false;
result.crossingSignature = zeros(1, 0);
result.numberOfPlaneCrossings = 0;
result.numberOfAcceptedCrossings = 0;

currentState = sourceState(:);
rejectedSinceAccepted = 0;

eventOptions = odeset(params.odeOptions, ...
    'Events', @(t, X) poincare_event(t, X, params));
noEventOptions = odeset(params.odeOptions, ...
    'Events', [], 'MaxStep', params.maxStepOffSection);

while result.numberOfAcceptedCrossings < horizon
    currentTime = 0;

    if abs(currentState(3) - params.zSection) <= params.sectionTolerance
        stepEnd = min( ...
            params.eventDisableTime, params.maxTimeToNextValidSection);
        if stepEnd <= 0
            return;
        end

        [~, XStep] = ode45( ...
            @(t, X) lorenz_rhs(t, X, params, pControl), ...
            [0, stepEnd], currentState, noEventOptions);
        currentState = XStep(end, :).';
        currentTime = stepEnd;
    end

    if currentTime >= params.maxTimeToNextValidSection
        return;
    end

    [~, ~, te, Xe] = ode45( ...
        @(t, X) lorenz_rhs(t, X, params, pControl), ...
        [currentTime, params.maxTimeToNextValidSection], ...
        currentState, eventOptions);

    if isempty(te)
        return;
    end

    currentState = Xe(end, :).';
    result.numberOfPlaneCrossings = result.numberOfPlaneCrossings + 1;

    if currentState(1) < params.xSectionMin
        result.crossingSignature(end + 1) = -1;
        rejectedSinceAccepted = rejectedSinceAccepted + 1;
        if rejectedSinceAccepted > params.maxRejectedCrossings
            return;
        end
    elseif currentState(1) > params.xSectionMin
        result.crossingSignature(end + 1) = +1;
        result.numberOfAcceptedCrossings = ...
            result.numberOfAcceptedCrossings + 1;
        rejectedSinceAccepted = 0;
    else
        % The branch boundary itself is neither an L nor an R crossing.
        result.crossingSignature(end + 1) = 0;
        rejectedSinceAccepted = rejectedSinceAccepted + 1;
        if rejectedSinceAccepted > params.maxRejectedCrossings
            return;
        end
    end
end

result.xFinal = currentState(1);
result.success = isfinite(result.xFinal);
end