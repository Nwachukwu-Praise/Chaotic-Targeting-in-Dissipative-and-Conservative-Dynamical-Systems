function results = test_shinbrot_replay_negative_control( ...
    sourceState, target, params, control)
%TEST_SHINBROT_REPLAY_NEGATIVE_CONTROL Prove the replay check can fail.
%
%   results = TEST_SHINBROT_REPLAY_NEGATIVE_CONTROL(sourceState, target, params, control)
%
%   Subsection 3.1.2 argues that replay verification is "a genuine check
%   rather than a restatement".  A 60/60 replay pass rate does not establish
%   that, because a check with no discriminating power also passes every
%   time.  This test supplies the missing evidence in two parts.
%
%   Part 1, independence.  For each search method it records the divergence
%   between the error the search reports and the error the replay measures.
%   A divergence of exactly zero means the two numbers were produced by the
%   same propagation from the same arguments, in which case the replay
%   cannot fail and verifies nothing.  This is presently the case for
%   hybridGridBisection: search_parameter_to_target sets finalTargetError
%   from x_after_n_crossings, and replay_selected_parameter re-runs the same
%   deterministic loop.  The test reports the fact rather than hiding it.
%
%   Part 2, negative control.  The stored perturbation and the stored
%   horizon are each deliberately corrupted and the replay is re-run.  A
%   replay with discriminating power must reject the corrupted candidate.
%   If a corrupted replay still passes, the check is not measuring what it
%   claims to measure and the test fails.
%
%   No production code path is modified.  The replay used here is the same
%   propagation the benchmark uses, so what is demonstrated applies directly
%   to the reported results.

methods = {verified_bisection_identifier(), 'hybridGridBisection'};
% A perturbation large enough to move the image well outside the target but
% small enough to stay inside the admissible interval, so the corrupted case
% remains a legitimate candidate rather than an out-of-bounds one.
perturbationOffsets = [1e-3, 5e-3];

rows = {};
fprintf('\n=== Replay independence and negative control ===\n');
fprintf('source X_s = %.8g,  eps_t = %.4g\n', sourceState(1), target.tolerance);

for m = 1:numel(methods)
    methodName = methods{m};
    methodControl = control;
    methodControl.searchMethod = methodName;

    search = search_parameter_to_target( ...
        sourceState, target, params, methodControl, methodName);

    if ~search.found
        fprintf('\n%s: no candidate found, skipped.\n', methodName);
        continue;
    end

    baseline = local_replay(sourceState, target, params, ...
        search.selectedP, search.horizon);
    divergence = abs(search.finalTargetError - baseline.targetError);
    independent = divergence > 0;

    fprintf('\n%s\n', methodName);
    fprintf('  selected p                     : %.10g\n', search.selectedP);
    fprintf('  horizon                        : %d\n', search.horizon);
    fprintf('  search-reported error          : %.12g\n', search.finalTargetError);
    fprintf('  replay-measured error          : %.12g\n', baseline.targetError);
    fprintf('  divergence                     : %.6g\n', divergence);
    fprintf('  baseline replay hit            : %d\n', baseline.hit);
    if independent
        fprintf('  independence                   : divergent, replay re-integrates\n');
    else
        fprintf(['  independence                   : ZERO divergence.  The search\n', ...
            '                                   error and the replay error come\n', ...
            '                                   from the same propagation, so the\n', ...
            '                                   agreement is guaranteed and is not\n', ...
            '                                   evidence of correctness.\n']);
    end

    rows(end + 1, :) = {string(methodName), "baseline", NaN, search.selectedP, ...
        search.horizon, baseline.targetError, baseline.hit, true, ...
        divergence, independent}; %#ok<AGROW>

    % ---- negative control 1: corrupt the stored perturbation ----------
    for k = 1:numel(perturbationOffsets)
        offset = perturbationOffsets(k);
        corruptedP = search.selectedP + offset;
        if abs(corruptedP) > control.deltaP
            corruptedP = search.selectedP - offset;
        end
        if abs(corruptedP) > control.deltaP
            continue;
        end
        corrupted = local_replay(sourceState, target, params, ...
            corruptedP, search.horizon);
        rejected = ~corrupted.hit;
        fprintf('  corrupt p by %+.4g            : replay hit = %d, error = %.6g -> %s\n', ...
            corruptedP - search.selectedP, corrupted.hit, corrupted.targetError, ...
            local_verdict(rejected));
        rows(end + 1, :) = {string(methodName), "corruptedParameter", offset, ...
            corruptedP, search.horizon, corrupted.targetError, corrupted.hit, ...
            rejected, NaN, independent}; %#ok<AGROW>
    end

    % ---- negative control 2: corrupt the stored horizon ---------------
    for shift = [1, -1]
        corruptedHorizon = search.horizon + shift;
        if corruptedHorizon < 1
            continue;
        end
        corrupted = local_replay(sourceState, target, params, ...
            search.selectedP, corruptedHorizon);
        rejected = ~corrupted.hit;
        fprintf('  corrupt horizon %+d            : replay hit = %d, error = %.6g -> %s\n', ...
            shift, corrupted.hit, corrupted.targetError, local_verdict(rejected));
        rows(end + 1, :) = {string(methodName), "corruptedHorizon", shift, ...
            search.selectedP, corruptedHorizon, corrupted.targetError, ...
            corrupted.hit, rejected, NaN, independent}; %#ok<AGROW>
    end
end

results = cell2table(rows, 'VariableNames', {'method','case','corruption', ...
    'replayedP','replayedHorizon','replayTargetError','replayHit', ...
    'expectedOutcomeObserved','searchReplayDivergence','searchIsIndependent'});

baselineRows = results.case == "baseline";
negativeRows = ~baselineRows;

allBaselinesPass = all(results.replayHit(baselineRows));
allNegativesRejected = all(results.expectedOutcomeObserved(negativeRows));

fprintf('\n--- summary ---\n');
disp(results);
fprintf('baseline replays passed            : %d\n', allBaselinesPass);
fprintf('corrupted replays rejected         : %d of %d\n', ...
    sum(results.expectedOutcomeObserved(negativeRows)), sum(negativeRows));
fprintf('methods with independent replay    : %d of %d\n', ...
    sum(results.searchIsIndependent(baselineRows)), sum(baselineRows));

if ~allBaselinesPass
    error('SchroerOtt:ReplayBaselineFailed', ...
        'A baseline replay failed; the reported candidate does not reproduce.');
end
if ~allNegativesRejected
    error('SchroerOtt:ReplayNegativeControlFailed', ...
        ['A deliberately corrupted candidate passed replay.  The replay ', ...
        'check has no discriminating power and must not be reported as ', ...
        'verification.']);
end
fprintf(['\nBoth negative controls behaved as required: the replay rejects\n', ...
    'a corrupted perturbation and a corrupted horizon.  Note separately\n', ...
    'that a zero search-replay divergence still means the reported error\n', ...
    'was not independently recomputed for that method.\n']);
fprintf('================================================\n');
end

% -------------------------------------------------------------------------

function replay = local_replay(sourceState, target, params, pControl, horizon)
%LOCAL_REPLAY Independent propagation, matching replay_selected_parameter.
replay.hit = false;
replay.targetError = Inf;
replay.failed = false;
replay.numCrossings = 0;

if ~isfinite(pControl) || ~isfinite(horizon) || horizon < 1
    replay.failed = true;
    return;
end

currentState = sourceState(:);
for crossing = 1:horizon
    segment = next_valid_section_crossing(currentState, params, pControl);
    if ~segment.success
        replay.failed = true;
        return;
    end
    currentState = segment.eventState;
    replay.numCrossings = crossing;
end

replay.targetError = abs(currentState(1) - target.x);
replay.hit = replay.targetError <= target.tolerance;
end

function s = local_verdict(rejected)
if rejected
    s = 'rejected as required';
else
    s = 'ACCEPTED, negative control FAILED';
end
end
