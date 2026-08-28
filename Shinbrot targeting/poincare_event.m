function [value, isterminal, direction] = poincare_event(~, X, params)
%POINCARE_EVENT Detect a crossing of the Poincare-section plane.
%
% The additional half-plane condition x > params.xSectionMin cannot be
% represented as a scalar zero event. next_valid_section_crossing therefore
% rejects plane crossings that do not satisfy it and continues integrating.

value = X(3) - params.zSection;
isterminal = 1;
direction = params.crossingDirection;
end
