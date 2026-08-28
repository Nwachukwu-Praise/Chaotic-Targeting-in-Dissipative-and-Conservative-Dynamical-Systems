function tf = so_connection_better(a, b)
%SO_CONNECTION_BETTER Deterministic connection tie-break.
if isempty(b) || ~isfield(b, 'success') || ~b.success
    tf = a.success;
    return;
end
if ~a.success
    tf = false;
    return;
end
keysA = [a.totalIterations, abs(a.control), a.intersectionResidual, ...
    phase_key(a.targetPhasePointIndex), chain_key(a.targetChainID), a.nForward];
keysB = [b.totalIterations, abs(b.control), b.intersectionResidual, ...
    phase_key(b.targetPhasePointIndex), chain_key(b.targetChainID), b.nForward];
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

