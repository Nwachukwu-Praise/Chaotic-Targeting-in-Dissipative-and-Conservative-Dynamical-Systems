function results = test_so_paper_ensemble_double_runtime()
%TEST_SO_PAPER_ENSEMBLE_DOUBLE_RUNTIME Focused checks for the versioned run.

audit = so_run_paper_ensemble_double_runtime('estimate', 'SaveOutputs', false, 'Quiet', true);
originalSignatures = original_file_signatures(audit.originalOutputFiles);

tests = {
    'new budget is exactly twice original', @test_budget_is_doubled
    'saved 50 source coordinates are reused exactly', @test_source_states_exact
    'scientific settings other than runtime are unchanged', @test_scientific_settings
    'short budget reports time budget exhausted', @test_short_budget_timeout
    'timeout is not numerical failure', @test_timeout_not_numerical_failure
    'checkpoint resume skips completed source', @test_checkpoint_resume
    'serial and parallel deterministic subset agree', @test_serial_parallel_subset
    'original ensemble files remain unchanged', @test_original_files_unchanged
    };

rows = {};
for i = 1:size(tests, 1)
    name = tests{i, 1};
    fn = tests{i, 2};
    timer = tic;
    skipped = false;
    message = "";
    try
        detail = fn();
        if isstruct(detail) && isfield(detail, 'skipped')
            skipped = detail.skipped;
        end
        passed = true;
    catch err
        passed = false;
        message = string(err.message);
    end
    rows(end + 1, :) = {string(name), passed, skipped, toc(timer), message}; %#ok<AGROW>
end

results = cell2table(rows, 'VariableNames', ...
    {'testName','passed','skipped','runtimeSeconds','message'});
disp(results);
if any(~results.passed)
    error('SchroerOtt:DoubleRuntimeTestsFailed', ...
        '%d doubled-runtime focused tests failed.', sum(~results.passed));
end

    function detail = test_budget_is_doubled()
        assert(audit.originalRuntimeBudgetSeconds == 120);
        assert(audit.newRuntimeBudgetSeconds == 2 * audit.originalRuntimeBudgetSeconds);
        assert(audit.newRuntimeBudgetSeconds == 240);
        detail.skipped = false;
    end

    function detail = test_source_states_exact()
        assert(audit.sourceStatesMatchOriginal);
        assert(isequaln(audit.sourceStates, audit.regeneratedSourceStates));
        assert(size(audit.sourceStates, 2) == 50);
        detail.skipped = false;
    end

    function detail = test_scientific_settings()
        assert(audit.scientificSettingsPreserved);
        assert(audit.newControlAmplitude == audit.originalControlAmplitude);
        assert(audit.independentLimits.maxControl == audit.originalControlAmplitude);
        assert(audit.independentLimits.maxForwardIterations == 10);
        assert(audit.independentLimits.maxBackwardIterations == 8);
        assert(audit.independentLimits.maxTotalTransferTime == 18);
        assert(isequal(audit.independentLimits.routeRotationNumbers, ...
            [1/4 1/3 2/5 1/2 3/5 2/3 3/4]));
        detail.skipped = false;
    end

    function detail = test_short_budget_timeout()
        outDir = fresh_output_dir('short_budget');
        pilot = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', 1, 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', outDir, 'OutputStem', 'unit_short_budget', ...
            'UseParallel', false, 'Quiet', true);
        assert(height(pilot.trials) == 1);
        assert(pilot.trials.terminationReason(1) == "time_budget_exhausted");
        detail.skipped = false;
    end

    function detail = test_timeout_not_numerical_failure()
        outDir = fresh_output_dir('timeout_class');
        pilot = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', 1, 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', outDir, 'OutputStem', 'unit_timeout_class', ...
            'UseParallel', false, 'Quiet', true);
        assert(pilot.trials.terminationReason(1) == "time_budget_exhausted");
        assert(pilot.trials.terminationReason(1) ~= "numerical_failure");
        detail.skipped = false;
    end

    function detail = test_checkpoint_resume()
        outDir = fresh_output_dir('checkpoint_resume');
        stem = 'unit_checkpoint_resume';
        first = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', 1, 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', outDir, 'OutputStem', stem, ...
            'UseParallel', false, 'Quiet', true);
        second = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', [1 2], 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', outDir, 'OutputStem', stem, ...
            'UseParallel', false, 'Quiet', true);
        assert(isequal(first.newlyCompletedIndices, 1));
        assert(isequal(second.skippedCompletedIndices, 1));
        assert(isequal(second.newlyCompletedIndices, 2));
        assert(height(second.trials) == 2);
        detail.skipped = false;
    end

    function detail = test_serial_parallel_subset()
        if audit.availableWorkers < 2
            detail.skipped = true;
            return;
        end
        serialDir = fresh_output_dir('serial_subset');
        parallelDir = fresh_output_dir('parallel_subset');
        serialRun = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', [1 2], 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', serialDir, 'OutputStem', 'unit_serial_subset', ...
            'UseParallel', false, 'Quiet', true);
        parallelRun = so_run_paper_ensemble_double_runtime('resumeDoubleRuntime', ...
            'SourceIndices', [1 2], 'RuntimeBudgetSeconds', 1e-6, ...
            'OutputDirectory', parallelDir, 'OutputStem', 'unit_parallel_subset', ...
            'UseParallel', true, 'NumWorkers', min(2, audit.availableWorkers), ...
            'Quiet', true);
        serialRows = sortrows(serialRun.trials(:, {'sourceIndex','sourceState', ...
            'found','terminationReason','iterations','controls','selectedRoute'}), 'sourceIndex');
        parallelRows = sortrows(parallelRun.trials(:, {'sourceIndex','sourceState', ...
            'found','terminationReason','iterations','controls','selectedRoute'}), 'sourceIndex');
        assert(isequaln(serialRows, parallelRows));
        detail.skipped = false;
    end

    function detail = test_original_files_unchanged()
        currentSignatures = original_file_signatures(audit.originalOutputFiles);
        assert(isequaln(originalSignatures, currentSignatures));
        detail.skipped = false;
    end
end

function folder = fresh_output_dir(label)
folder = [tempname, '_', label];
mkdir(folder);
end

function signatures = original_file_signatures(files)
files = string(files(:));
bytes = zeros(numel(files), 1);
datenumValue = zeros(numel(files), 1);
for i = 1:numel(files)
    info = dir(files(i));
    if isempty(info)
        error('SchroerOtt:MissingOriginalOutput', ...
            'Original output file is missing: %s', files(i));
    end
    bytes(i) = info.bytes;
    datenumValue(i) = info.datenum;
end
signatures = table(files, bytes, datenumValue, ...
    'VariableNames', {'file','bytes','datenumValue'});
end
