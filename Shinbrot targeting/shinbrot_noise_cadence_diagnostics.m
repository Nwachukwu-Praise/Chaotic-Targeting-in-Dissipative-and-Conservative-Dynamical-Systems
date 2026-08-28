function diagnostics = shinbrot_noise_cadence_diagnostics(mapData, noise, burnIn)
%SHINBROT_NOISE_CADENCE_DIAGNOSTICS Measure retarget cadence against map time.

if nargin < 3 || isempty(burnIn)
    burnIn = 100;
end
if ~isfield(mapData, 'eventTimes') || numel(mapData.eventTimes) < burnIn + 2
    error('mapData.eventTimes does not contain enough return times.');
end

eventTimes = mapData.eventTimes(:);
intervals = diff(eventTimes((burnIn + 1):end));
intervals = intervals(isfinite(intervals) & intervals > 0);
if isempty(intervals)
    error('No positive finite return intervals were available.');
end

diagnostics.burnIn = burnIn;
diagnostics.Tbar = mean(intervals);
diagnostics.medianReturnTime = median(intervals);
diagnostics.dt = noise.dt;
diagnostics.stepsPerMeanReturn = diagnostics.Tbar / noise.dt;
diagnostics.retargetStepInterval = noise.retargetStepInterval;
diagnostics.retargetPhysicalTime = ...
    noise.retargetStepInterval * noise.dt;
diagnostics.cadenceRatio = ...
    diagnostics.retargetPhysicalTime / diagnostics.Tbar;
diagnostics.paperReferenceStepCount = ...
    get_field_with_default(noise, 'paperReferenceStepCount', NaN);
diagnostics.retargetCadencePolicy = ...
    string(get_field_with_default(noise, 'retargetCadencePolicy', ''));
diagnostics.literalFortyStepPhysicalTime = 40 * noise.dt;
diagnostics.literalFortyStepCadenceRatio = ...
    diagnostics.literalFortyStepPhysicalTime / diagnostics.Tbar;

if isfield(noise, 'maxControlledCrossings') && ...
        isfield(noise, 'maxRk4Steps')
    requiredSteps = noise.maxControlledCrossings * ...
        diagnostics.stepsPerMeanReturn;
    diagnostics.configuredSafetyFactor = ...
        noise.maxRk4Steps / requiredSteps;
else
    diagnostics.configuredSafetyFactor = NaN;
end

tableStruct = diagnostics;
diagnostics.table = struct2table(tableStruct, 'AsArray', true);
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
