function times = so_uncontrolled_transport_time(sourceStates, cfg, iterationCap)
%SO_UNCONTROLLED_TRANSPORT_TIME Natural first-passage time into the target.
%
%   times = SO_UNCONTROLLED_TRANSPORT_TIME(Z0, cfg, cap)
%
% Z0 is 2-by-N.  Returns, for each column, the number of iterations of the
% UNCONTROLLED standard map before the orbit first lands in
% cfg.targetRectangle, or NaN if it has not done so within cap iterations.
%
% This is the comparison that gives the controlled numbers their meaning.
% Schroer and Ott quote 1119 to 3.77e6 uncontrolled steps against 125 to 132
% controlled ones for the same source region; without the uncontrolled
% figure, a controlled transfer time is just a number.
%
% All N orbits are advanced together, so the cost is cap vector steps rather
% than N * cap scalar steps.  Orbits that have already arrived are frozen so
% a long-lived one cannot corrupt an early arrival.
if nargin < 3 || isempty(iterationCap)
    iterationCap = 4e6;
end
rect = cfg.targetRectangle;
rect.yMin = rect.yMin + cfg.transport.targetLiftShift;
rect.yMax = rect.yMax + cfg.transport.targetLiftShift;

n = size(sourceStates, 2);
Z = sourceStates;
times = nan(1, n);
active = true(1, n);

% A point already sitting in the target has transport time zero.
hit = so_point_in_rectangle(Z, rect, 0);
times(hit) = 0;
active(hit) = false;

blockSize = 20000;
done = 0;
while done < iterationCap && any(active)
    thisBlock = min(blockSize, iterationCap - done);
    idx = find(active);
    Zi = Z(:, idx);
    for step = 1:thisBlock
        Zi = so_standard_map_lifted(Zi, cfg);
        inside = so_point_in_rectangle(Zi, rect, 0);
        if any(inside)
            arrived = idx(inside);
            times(arrived) = done + step;
            active(arrived) = false;
            keep = ~inside;
            Zi = Zi(:, keep);
            idx = idx(keep);
            if isempty(idx)
                break;
            end
        end
    end
    Z(:, idx) = Zi;
    done = done + thisBlock;
end
end
