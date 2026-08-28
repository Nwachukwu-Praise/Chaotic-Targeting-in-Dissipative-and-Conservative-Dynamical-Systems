function parallel = select_noise_parallel_workers(requestedWorkers)
%SELECT_NOISE_PARALLEL_WORKERS Choose a safe process-pool size.
%
% requestedWorkers may be 'auto' or a positive integer.  The selected count
% never exceeds the local process-worker limit.  GPU hardware is deliberately
% ignored because the noisy Lorenz workflow is branch-heavy, three-state RK4
% code with bisection retargeting rather than a dense array calculation.

if nargin < 1 || isempty(requestedWorkers)
    requestedWorkers = 'auto';
end

parallel = struct();
parallel.requestedWorkers = requestedWorkers;
parallel.hasParallelToolbox = license('test', 'Distrib_Computing_Toolbox');
parallel.availableWorkers = 1;
parallel.existingPoolWorkers = 0;
parallel.existingPoolClass = '';
parallel.candidateWorkers = 1;
parallel.selectedWorkers = 1;
parallel.executionMode = 'serial';
parallel.gpuAccelerationExcluded = true;
parallel.selectionReason = 'Parallel Computing Toolbox unavailable.';
parallel.parallelEfficiencyAssumption = 0.70;

if ~parallel.hasParallelToolbox
    return;
end

try
    localCluster = parcluster('local');
    parallel.availableWorkers = max(1, localCluster.NumWorkers);
catch problem
    parallel.selectionReason = ['Could not inspect the local cluster: ', ...
        problem.message];
    parallel.hasParallelToolbox = false;
    parallel.availableWorkers = 1;
    return;
end

try
    pool = gcp('nocreate');
    if ~isempty(pool)
        parallel.existingPoolWorkers = pool.NumWorkers;
        parallel.existingPoolClass = class(pool);
    end
catch
    parallel.existingPoolWorkers = 0;
    parallel.existingPoolClass = '';
end

candidates = unique([4, 8, 12, 16, parallel.availableWorkers]);
candidates = candidates(candidates >= 1 & ...
    candidates <= parallel.availableWorkers);
if isempty(candidates)
    candidates = 1;
end
parallel.candidateWorkers = candidates;

if isnumeric(requestedWorkers)
    requested = round(requestedWorkers);
    if ~isscalar(requested) || ~isfinite(requested) || requested < 1
        error('parallel.requestedWorkers must be ''auto'' or a positive integer.');
    end
    parallel.selectedWorkers = min(requested, parallel.availableWorkers);
    if requested > parallel.availableWorkers
        parallel.selectionReason = sprintf( ...
            'Requested %d workers but the local process-worker limit is %d.', ...
            requested, parallel.availableWorkers);
    else
        parallel.selectionReason = sprintf( ...
            'Using explicitly requested process-worker count: %d.', ...
            parallel.selectedWorkers);
    end
else
    requestedText = char(string(requestedWorkers));
    if ~strcmpi(requestedText, 'auto')
        error('parallel.requestedWorkers must be ''auto'' or a positive integer.');
    end

    if parallel.availableWorkers > 2
        parallel.selectedWorkers = parallel.availableWorkers - 1;
    else
        parallel.selectedWorkers = parallel.availableWorkers;
    end
    parallel.selectionReason = sprintf([ ...
        'Auto-selected close to the process-worker limit while leaving ', ...
        'capacity for Windows/MATLAB when possible.']);
end

if parallel.selectedWorkers > 1
    parallel.executionMode = 'parallel';
else
    parallel.executionMode = 'serial';
end
end
