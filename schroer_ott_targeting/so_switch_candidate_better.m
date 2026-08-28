function tf = so_switch_candidate_better(a, b)
%SO_SWITCH_CANDIDATE_BETTER Tie-break for J(j) probes.
if isempty(b) || ~isfield(b, 'finite') || ~b.finite
    tf = a.finite;
    return;
end
if ~a.finite
    tf = false;
    return;
end
keysA = [a.objectiveJ, a.j, abs(a.connection.control), ...
    a.connection.intersectionResidual, phase_key(a.connection.targetPhasePointIndex), ...
    chain_key(a.connection.targetChainID)];
keysB = [b.objectiveJ, b.j, abs(b.connection.control), ...
    b.connection.intersectionResidual, phase_key(b.connection.targetPhasePointIndex), ...
    chain_key(b.connection.targetChainID)];
tf = false;
for i = 1:numel(keysA)
    if keysA(i) < keysB(i) - 1e-14
        tf = true;
        return;
    elseif keysA(i) > keysB(i) + 1e-14
        return;
    end
end
end

function v = phase_key(x)
if isnan(x)
    v = Inf;
else
    v = x;
end
end

function v = chain_key(id)
txt = char(id);
nums = regexp(txt, '\d+', 'match');
if isempty(nums)
    v = Inf;
else
    v = str2double(nums{end});
end
end

