function config = build_bisection_noise_config( ...
    sourceState, target, params, noiseControl, noise, sweep)
%BUILD_BISECTION_NOISE_CONFIG Scientific configuration for noise results.

verifiedIdentifier = verified_bisection_identifier();
if ~isfield(noiseControl, 'searchMethod') || ...
        ~strcmp(noiseControl.searchMethod, verifiedIdentifier)
    error('noiseControl.searchMethod must be %s.', verifiedIdentifier);
end
if ~isfield(sweep, 'searchMethod') || ...
        ~strcmp(sweep.searchMethod, verifiedIdentifier)
    error('sweep.searchMethod must be %s.', verifiedIdentifier);
end
if ~isfield(sweep, 'modes') || ~all_valid_modes(sweep.modes)
    error('The noise sweep defines an unknown retargeting mode.');
end
if ~isequal(sweep.modes(:).', {'paperCycleCadenceRetargetMode'})
    error(['The final repaired noise experiment must use only ', ...
        'paperCycleCadenceRetargetMode.']);
end
if ~isfield(sweep, 'sigmaNoiseValues')
    error('The final noise sweep must define sweep.sigmaNoiseValues.');
end
expectedSigmaNoiseValues = [0, 0.01, 0.03, 0.05, 0.075, 0.1, ...
    0.15, 0.2, 0.3, 0.5, 1.0];
if ~isequal(sweep.sigmaNoiseValues(:).', expectedSigmaNoiseValues)
    error(['The final noise sweep must use sigmaNoiseValues = ', ...
        mat2str(expectedSigmaNoiseValues), '.']);
end
if ~isfield(sweep, 'numTrials') || sweep.numTrials ~= 10
    error('The final noise sweep must use sweep.numTrials = 10.');
end
if ~isfield(noise, 'retargetCadencePolicy') || ...
        ~strcmp(noise.retargetCadencePolicy, 'meanReturnCycleMatched')
    error(['The repaired final noise workflow requires ', ...
        'noise.retargetCadencePolicy = ''meanReturnCycleMatched''.']);
end
if ~isfield(noise, 'paperReferenceStepCount') || ...
        noise.paperReferenceStepCount ~= 40
    error('The repaired final noise workflow requires paperReferenceStepCount = 40.');
end
if ~isfield(noise, 'retargetStepInterval') || ...
        ~isfinite(noise.retargetStepInterval) || ...
        noise.retargetStepInterval < 1
    error('The repaired final noise workflow requires a finite retargetStepInterval.');
end
if ~isfield(noise, 'cadenceRatio') || abs(noise.cadenceRatio - 1) > 0.05
    error('The repaired final noise workflow requires abs(cadenceRatio - 1) <= 0.05.');
end
if ~isfield(noise, 'maxRk4Steps') || ~isfinite(noise.maxRk4Steps)
    error('noise.maxRk4Steps must be a finite global full-trial budget.');
end
if ~isfield(noise, 'maxControlledCrossings') || ...
        noise.maxControlledCrossings ~= 60
    error('The final noise workflow requires noise.maxControlledCrossings = 60.');
end
if isfield(noise, 'requiredRk4Steps') && ...
        noise.maxRk4Steps < noise.requiredRk4Steps
    error('noise.maxRk4Steps is below the measured required RK4 budget.');
end
if ~isfield(sweep, 'isSmokeTest') || sweep.isSmokeTest
    error('Full noise configuration requires sweep.isSmokeTest == false.');
end

config.verifiedSearchMethod = verifiedIdentifier;
% The dispatcher identifier names the case, not the function behind it, so
% it does not change when the search implementation changes.  The
% fingerprint below does, which means saved results produced by a different
% search are rejected instead of silently reused.
config.searchImplementation = verified_bisection_implementation();
config.noiseImplementation = shinbrot_dependency_fingerprint('noise');
config.allImplementation = shinbrot_dependency_fingerprint('all');
config.modes = sweep.modes(:).';
config.retargetingMode = strjoin(config.modes, ', ');
config.sigmaNoiseValues = sweep.sigmaNoiseValues(:).';
config.numTrials = sweep.numTrials;
config.baseSeed = sweep.baseSeed;
config.seedFormula = 'baseSeed + 100000*modeIndex + 1000*sigmaIndex + trialIndex';
if isfield(sweep, 'sourceIndex')
    config.sourceIndex = sweep.sourceIndex;
else
    config.sourceIndex = NaN;
end
config.sourceState = sourceState(:);
config.targetState = target.state(:);
config.targetX = target.x;
config.targetTolerance = target.tolerance;
config.deltaP = noiseControl.deltaP;
config.maxSearchCrossings = noiseControl.maxSearchCrossings;
config.bisectionIterations = noiseControl.bisectionIterations;
config.maxDiscontinuityIsolationDepth = ...
    get_field_with_default(noiseControl, 'maxDiscontinuityIsolationDepth', NaN);
config.maxDiscontinuityIntervals = ...
    get_field_with_default(noiseControl, 'maxDiscontinuityIntervals', NaN);
config.dt = noise.dt;
config.maxRk4Steps = get_field_with_default(noise, 'maxRk4Steps', NaN);
config.meanReturnTimeForBudget = ...
    get_field_with_default(noise, 'meanReturnTimeForBudget', NaN);
config.estimatedStepsPerCrossing = ...
    get_field_with_default(noise, 'estimatedStepsPerCrossing', NaN);
config.budgetSafetyFactor = ...
    get_field_with_default(noise, 'budgetSafetyFactor', NaN);
config.requiredRk4Steps = ...
    get_field_with_default(noise, 'requiredRk4Steps', NaN);
config.budgetBurnIn = get_field_with_default(noise, 'budgetBurnIn', NaN);
config.paperReferenceStepCount = ...
    get_field_with_default(noise, 'paperReferenceStepCount', NaN);
config.retargetCadencePolicy = ...
    get_field_with_default(noise, 'retargetCadencePolicy', '');
config.retargetStepInterval = ...
    get_field_with_default(noise, 'retargetStepInterval', NaN);
config.retargetPhysicalTime = ...
    get_field_with_default(noise, 'retargetPhysicalTime', NaN);
config.cadenceRatio = get_field_with_default(noise, 'cadenceRatio', NaN);
config.maxTimeToNextValidSection = noise.maxTimeToNextValidSection;
config.maxControlledCrossings = noise.maxControlledCrossings;
config.maxRejectedCrossings = noise.maxRejectedCrossings;
config.maxStateNorm = noise.maxStateNorm;
config.zSection = params.zSection;
config.xSectionMin = params.xSectionMin;
config.crossingDirection = params.crossingDirection;
config.sectionTolerance = params.sectionTolerance;
config.sigma = params.sigma;
config.rho = params.rho;
config.beta = params.beta;
config.noiseUpdate = ...
    'stepStartState = currentState + sigmaNoise*randn(3,1); stepEndState = rk4_lorenz_step(stepStartState,dt,params,pControl)';
config.crossingInterpolation = ...
    'linear interpolation between stepStartState and stepEndState; instantaneous noise jumps excluded';
config.totalTrials = numel(config.sigmaNoiseValues) * ...
    config.numTrials * numel(sweep.modes);
if isfield(sweep, 'totalTrials') && sweep.totalTrials ~= config.totalTrials
    error('sweep.totalTrials does not match the executable configuration.');
end
config.isSmokeTest = false;
config.trajectoryRetentionPolicy = get_field_with_default(sweep, ...
    'trajectoryRetentionPolicy', 'not declared');
config.signatureText = stable_config_text(config);
end

function tf = all_valid_modes(modes)
validModes = {'paperCycleCadenceRetargetMode', ...
    'paper40StepRetargetMode', 'sectionRetargetMode'};
if ~iscell(modes)
    tf = false;
    return;
end
tf = all(ismember(modes(:), validModes));
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end

function textValue = stable_config_text(config)
fields = sort(fieldnames(config));
parts = cell(numel(fields), 1);
for i = 1:numel(fields)
    name = fields{i};
    if strcmp(name, 'signatureText')
        continue;
    end
    value = config.(name);
    if isnumeric(value) || islogical(value)
        valueText = mat2str(value, 17);
    elseif ischar(value)
        valueText = value;
    elseif isstring(value)
        valueText = char(value);
    else
        valueText = evalc('disp(value)');
    end
    parts{i} = [name, '=', valueText];
end
textValue = strjoin(parts(~cellfun('isempty', parts)), newline);
end
