function convergence = fixed_step_open_loop_convergence( ...
    sourceState, target, params, noise, pControl, horizon, dtValues)
%FIXED_STEP_OPEN_LOOP_CONVERGENCE Open-loop zero-noise RK4 dt study.

if nargin < 7 || isempty(dtValues)
    dtValues = noise.dt ./ [1, 2, 4];
end
dtValues = dtValues(:);
rows = cell(numel(dtValues), 1);

for k = 1:numel(dtValues)
    trialNoise = noise;
    trialNoise.dt = dtValues(k);
    nominalTimeBudget = get_field_with_default(noise, 'maxRk4Steps', 1) * ...
        noise.dt;
    trialNoise.maxRk4Steps = max(ceil(nominalTimeBudget / dtValues(k)), ...
        ceil(20 / dtValues(k)));
    comparison = fixed_step_open_loop_replay(sourceState, target, params, ...
        trialNoise, pControl, horizon);
    rows{k} = comparison.table;
end

convergence = vertcat(rows{:});
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
