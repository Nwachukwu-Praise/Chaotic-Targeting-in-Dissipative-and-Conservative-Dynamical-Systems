function trial = initialize_noisy_trial_record( ...
    sourceState, target, sigmaNoise, trialNumber, seed, mode, ...
    verifiedIdentifier, isSmokeTest)
%INITIALIZE_NOISY_TRIAL_RECORD Shared schema for noisy targeting trials.

trial.sigmaNoise = sigmaNoise;
trial.trialNumber = trialNumber;
trial.seed = seed;
trial.mode = mode;
trial.searchMethod = verifiedIdentifier;
trial.isSmokeTest = isSmokeTest;

trial.hit = false;
trial.numCrossings = NaN;             % Backward-compatible alias.
trial.numCrossingsToTarget = NaN;
trial.numAcceptedCrossings = 0;
trial.numRejectedCrossings = 0;
trial.physicalTime = NaN;
trial.numRk4Steps = 0;
trial.integrationFailed = false;
trial.stepBudgetExhausted = false;
trial.acceptedCrossingLimitExhausted = false;
trial.rejectedCrossingLimitExhausted = false;
trial.failureReason = '';

trial.initialSearchPerformed = false;
trial.initialSearch = struct();
trial.initialSearchMethod = verifiedIdentifier;
trial.initialSearchSelectionMethod = '';
trial.initialSearchGridEvaluations = NaN;
trial.initialSearchCrossingOrderTestUsed = false;
trial.initialSearchFound = false;
trial.initialSearchSelectedP = NaN;
trial.initialSearchHorizon = NaN;
trial.initialSearchFinalTargetError = NaN;
trial.numSearchEvaluations = 0;
trial.numCompleteBisectionSearches = 0;
trial.numRetargetingEvents = 0;
trial.numRetargetingFailures = 0;

trial.selectedPInitial = NaN;
trial.selectedPFinal = NaN;
trial.selectedPValues = [];
trial.finalDistance = abs(sourceState(1) - target.x);
trial.finalState = sourceState(:);
trial.states = [];
trial.xSequence = sourceState(1);
trial.eventTimes = [];
trial.t = [];
trial.X = [];
trial.computationalTime = NaN;

trial.retargetingCrossingIndices = [];
trial.retargetingStepIndices = [];
trial.retargetingTimes = [];
trial.retargetingSearchMethods = {};
trial.retargetingSelectedP = [];
trial.retargetingSearchFound = [];
trial.retargetingSelectionMethods = {};
trial.retargetingGridEvaluations = [];
trial.retargetingCrossingOrderTestUsed = [];

trial.acceptedCrossingStates = [];
trial.acceptedCrossingTimes = [];
trial.acceptedCrossingStepIndices = [];
trial.acceptedCrossingX = [];
trial.rejectedCrossingStates = [];
trial.rejectedCrossingTimes = [];
trial.rejectedCrossingStepIndices = [];
trial.retargetStepInterval = NaN;
trial.retargetCadencePolicy = '';
trial.paperReferenceStepCount = NaN;
trial.retargetPhysicalTime = NaN;
trial.cadenceRatio = NaN;
end
