function trial = run_noisy_targeting_trial( ...
    sourceState, target, params, noiseControl, noise, ...
    sigmaNoise, trialNumber, seed, mode, initialSearch, isSmokeTest)
%RUN_NOISY_TARGETING_TRIAL Dispatch one bisection-only noisy trial.
%
% paperCycleCadenceRetargetMode is the final paper-mode reconstruction:
% Gaussian coordinate noise is added before each fixed RK4 step, and the
% targeting calculation is repeated after about one accepted return-map
% cycle.  paper40StepRetargetMode is retained only as a literal historical
% audit path.  sectionRetargetMode is a labelled project extension that
% retargets only after accepted Poincare crossings.

verifiedIdentifier = verified_bisection_identifier();
noiseControl = require_bisection_control(noiseControl, verifiedIdentifier);

if nargin < 9 || isempty(mode)
    mode = 'paperCycleCadenceRetargetMode';
end
mode = normalize_mode(mode);
if nargin < 10 || isempty(initialSearch)
    initialSearch = search_parameter_to_target( ...
        sourceState, target, params, noiseControl, verifiedIdentifier);
end
assert_verified_bisection_search( ...
    initialSearch, verifiedIdentifier, 'initial noise-trial search', target);
if nargin < 11 || isempty(isSmokeTest)
    isSmokeTest = false;
end

rng(seed);
trialTimer = tic;
trial = initialize_noisy_trial_record(sourceState, target, sigmaNoise, ...
    trialNumber, seed, mode, verifiedIdentifier, isSmokeTest);

initialAccepted = initialSearch.found && isfinite(initialSearch.selectedP);
trial = record_noisy_search_event( ...
    trial, 'initial', 0, 0, 0, initialSearch, initialAccepted);

if ~initialAccepted
    trial.integrationFailed = true;
    trial.failureReason = 'initial deterministic bisection search failed';
    trial.computationalTime = toc(trialTimer);
    return;
end

pCurrent = initialSearch.selectedP;
trial.selectedPInitial = pCurrent;
trial.selectedPFinal = pCurrent;
trial.selectedPValues = pCurrent;

trialNoise = noise;
trialNoise.sigmaNoise = sigmaNoise;
trialNoise.noiseCoordinateSigma = sigmaNoise;
trialNoise.delta = sigmaNoise;

switch mode
    case 'paperCycleCadenceRetargetMode'
        trial = run_noisy_trial_paper_cycle_cadence( ...
            trial, sourceState, target, params, noiseControl, ...
            trialNoise, pCurrent, verifiedIdentifier, trialTimer);
    case 'paper40StepRetargetMode'
        trial = run_noisy_trial_paper40_step( ...
            trial, sourceState, target, params, noiseControl, ...
            trialNoise, pCurrent, verifiedIdentifier, trialTimer);
    case 'sectionRetargetMode'
        trial = run_noisy_trial_section_retarget( ...
            trial, sourceState, target, params, noiseControl, ...
            trialNoise, pCurrent, verifiedIdentifier, trialTimer);
    otherwise
        error('Unknown noisy retargeting mode: %s', mode);
end
end

function noiseControl = require_bisection_control(noiseControl, verifiedIdentifier)
if ~isfield(noiseControl, 'searchMethod') || ...
        ~strcmp(noiseControl.searchMethod, verifiedIdentifier)
    error('noiseControl.searchMethod must be %s.', verifiedIdentifier);
end
end

function mode = normalize_mode(mode)
if isa(mode, 'string')
    mode = char(mode);
end
validModes = {'paperCycleCadenceRetargetMode', ...
    'paper40StepRetargetMode', 'sectionRetargetMode'};
if ~any(strcmp(mode, validModes))
    error('Unknown noisy retargeting mode: %s', mode);
end
end
