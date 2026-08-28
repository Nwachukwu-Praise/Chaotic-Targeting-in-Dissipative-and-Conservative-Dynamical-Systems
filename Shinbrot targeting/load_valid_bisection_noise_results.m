function [noiseResults, isValid, differences] = load_valid_bisection_noise_results( ...
    fileName, currentConfig)
%LOAD_VALID_BISECTION_NOISE_RESULTS Load only a complete bisection noise result.

noiseResults = struct();
isValid = false;
differences = strings(0, 1);
if ~isfile(fileName)
    differences(end + 1, 1) = "file not found";
    return;
end

saved = load(fileName);
if isfield(saved, 'noiseResults')
    candidate = saved.noiseResults;
elseif isfield(saved, 'results')
    candidate = saved.results;
else
    differences(end + 1, 1) = "no noiseResults variable";
    return;
end

requiredFields = {'trials', 'trialsTable', 'summaryTable', ...
    'configuration', 'completed', 'isSmokeTest', 'searchMethod', ...
    'sigmaNoiseValues', 'numTrials', 'modes', 'totalTrials'};
if ~isstruct(candidate) || ~all(isfield(candidate, requiredFields))
    differences(end + 1, 1) = "missing required result fields";
    return;
end

[configValid, configDifferences] = validate_bisection_noise_config( ...
    candidate.configuration, currentConfig);
if ~configValid
    differences = [differences; configDifferences(:)];
    return;
end
verifiedIdentifier = currentConfig.verifiedSearchMethod;
expectedTotal = currentConfig.totalTrials;

if candidate.isSmokeTest
    differences(end + 1, 1) = "result is a smoke test";
end
if isfield(candidate, 'isPilotOnly') && candidate.isPilotOnly
    differences(end + 1, 1) = "result is pilot-only";
end
if ~candidate.completed
    differences(end + 1, 1) = "result is incomplete";
end
if ~strcmp(candidate.searchMethod, verifiedIdentifier)
    differences(end + 1, 1) = "result search method is not verified bisection";
end
if ~isequal(candidate.modes(:), currentConfig.modes(:))
    differences(end + 1, 1) = "result retargeting modes do not match current configuration";
end
if candidate.totalTrials ~= expectedTotal || numel(candidate.trials) ~= expectedTotal
    differences(end + 1, 1) = "result trial count mismatch";
end
if expectedTotal ~= 110
    differences(end + 1, 1) = "current configuration is not the repaired 110-trial experiment";
end
if ~istable(candidate.summaryTable) || ~istable(candidate.trialsTable)
    differences(end + 1, 1) = "summary or trial output is not a table";
end
if istable(candidate.trialsTable) && ...
        any(strcmp(candidate.trialsTable.Properties.VariableNames, ...
        'SearchMethod')) && ...
        any(strcmp(candidate.trialsTable.SearchMethod, 'hybridGridBisection'))
    differences(end + 1, 1) = "trial table contains hybridGridBisection";
end

for i = 1:numel(candidate.trials)
    try
        validate_trial(candidate.trials(i), verifiedIdentifier, currentConfig.modes);
    catch problem
        differences(end + 1, 1) = "trial " + i + ": " + problem.message;
        break;
    end
end

if isfield(candidate, 'completed') && candidate.completed && ...
        isfield(candidate, 'trials')
    completedModes = {candidate.trials.mode};
    if ~all(strcmp(completedModes, 'paperCycleCadenceRetargetMode'))
        differences(end + 1, 1) = "not all completed trials use paperCycleCadenceRetargetMode";
    end
end

isValid = isempty(differences);
if isValid
    noiseResults = candidate;
end
end

function validate_trial(trial, verifiedIdentifier, expectedModes)
if ~isfield(trial, 'searchMethod') || ~strcmp(trial.searchMethod, verifiedIdentifier)
    error('search method mismatch');
end
if ~isfield(trial, 'mode') || ~any(strcmp(trial.mode, expectedModes))
    error('retargeting mode mismatch');
end
if ~isfield(trial, 'isSmokeTest') || trial.isSmokeTest
    error('smoke-test trial');
end
if ~isfield(trial, 'initialSearchPerformed') || ~trial.initialSearchPerformed
    error('missing initial search provenance');
end
if ~isfield(trial, 'initialSearchMethod') || ...
        ~strcmp(trial.initialSearchMethod, verifiedIdentifier)
    error('initial search method mismatch');
end
if isfield(trial, 'initialSearchGridEvaluations') && ...
        trial.initialSearchGridEvaluations ~= 0
    error('initial search contains grid evaluations');
end
if ~isfield(trial, 'retargetingSearchMethods')
    error('missing retargeting provenance field');
end
if any(strcmp(trial.retargetingSearchMethods, 'hybridGridBisection')) || ...
        (~isempty(trial.retargetingSearchMethods) && ...
        ~all(strcmp(trial.retargetingSearchMethods, verifiedIdentifier)))
    error('non-bisection retargeting provenance');
end
if isfield(trial, 'retargetingGridEvaluations') && ...
        any(trial.retargetingGridEvaluations ~= 0)
    error('hybrid grid evaluations present');
end
if isfield(trial, 'retargetingSelectionMethods') && ...
        any(contains(trial.retargetingSelectionMethods, 'grid-first'))
    error('hybrid grid-first selection present');
end
if strcmp(trial.mode, 'paperCycleCadenceRetargetMode')
    if ~isfield(trial, 'retargetCadencePolicy') || ...
            ~strcmp(trial.retargetCadencePolicy, 'meanReturnCycleMatched')
        error('cycle-cadence trial policy mismatch');
    end
    if ~isfield(trial, 'paperReferenceStepCount') || ...
            trial.paperReferenceStepCount ~= 40
        error('cycle-cadence trial does not preserve paper reference step count');
    end
    if ~isfield(trial, 'retargetStepInterval') || ...
            ~isfinite(trial.retargetStepInterval) || ...
            trial.retargetStepInterval < 1
        error('cycle-cadence trial has invalid retargetStepInterval');
    end
    if ~isfield(trial, 'cadenceRatio') || abs(trial.cadenceRatio - 1) > 0.05
        error('cycle-cadence trial cadence ratio mismatch');
    end
    if isfield(trial, 'numRetargetingEvents') && ...
            isfield(trial, 'numRk4Steps') && ...
            trial.numRetargetingEvents > floor(trial.numRk4Steps / ...
            trial.retargetStepInterval)
        error('cycle-cadence retarget count exceeds scheduled cadence');
    end
elseif strcmp(trial.mode, 'paper40StepRetargetMode')
    if ~isfield(trial, 'retargetStepInterval') || trial.retargetStepInterval ~= 40
        error('paper mode retargetStepInterval is not 40');
    end
    if isfield(trial, 'retargetingStepIndices') && ...
            ~isempty(trial.retargetingStepIndices) && ...
            (~all(mod(trial.retargetingStepIndices, 40) == 0) || ...
            (numel(trial.retargetingStepIndices) > 1 && ...
            any(diff(trial.retargetingStepIndices) ~= 40)))
        error('paper-mode retargeting cadence mismatch');
    end
end
end
