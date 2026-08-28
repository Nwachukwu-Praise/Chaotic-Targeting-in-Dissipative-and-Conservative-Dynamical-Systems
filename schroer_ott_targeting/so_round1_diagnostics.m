function so_round1_diagnostics(runHotspots, runStages)
%SO_ROUND1_DIAGNOSTICS Correction round 1: verify, then diagnose.
%
%   Run this from the schroer_ott_targeting folder:
%
%       so_round1_diagnostics                % census only, fastest
%       so_round1_diagnostics(true, false)   % + MATLAB code profiler
%       so_round1_diagnostics(true, true)    % + one full trial, 30 min cap
%
%   Everything printed is written to so_round1_diagnostics.log as well.
%
%   Part A verifies that the rewritten so_build_curve reproduces the
%   previous builder exactly, and reports the speed-up.  Nothing further is
%   trustworthy unless part A passes.
%
%   Part B checks the newly stored monodromy spectrum and draws the
%   eigendirection figure.
%
%   Part C measures where the paper-geometry trial spends its time.  The
%   initial 50-source ensemble returned 0/50 within its declared budget:
%   45 sources hit the 120 s cap, 3 ended at provisional connection search,
%   and 2 exhausted switch-probe searches.  The outcome is computationally
%   unresolved at that budget, not a mathematical failure claim.

if nargin < 1 || isempty(runHotspots)
    runHotspots = false;
end
if nargin < 2 || isempty(runStages)
    runStages = false;
end

here = fileparts(mfilename('fullpath'));
addpath(here);
cd(here);

logFile = fullfile(here, 'so_round1_diagnostics.log');
if isfile(logFile)
    delete(logFile);
end
diary(logFile);
diary on;
cleanup = onCleanup(@() diary('off'));

fprintf('===============================================================\n');
fprintf('Schroer-Ott targeting - correction round 1 diagnostics\n');
fprintf('date            : %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>
fprintf('MATLAB version  : %s\n', version);
fprintf('folder          : %s\n', here);
fprintf('Mapping Toolbox : %d\n', license('test', 'MAP_Toolbox'));
fprintf('===============================================================\n');

% ---------------------------------------------------------------- part A
fprintf('\n\n############ PART A: so_build_curve equivalence ############\n');
try
    equivalence = test_so_build_curve_equivalence(); %#ok<NASGU>
    fprintf('PART A: PASSED\n');
catch err
    fprintf('PART A: FAILED - %s\n', err.message);
    fprintf('%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
    diary off;
    return;
end

% ---------------------------------------------------------------- part B
fprintf('\n\n############ PART B: monodromy spectrum ############\n');
cfg = schroer_ott_default_config();
cfg.verbose = false;
catalogue = so_enumerate_periodic_orbits(cfg);
chains = catalogue.chains;

eigenTable = so_chain_eigen_table(chains, false);
fprintf('\nHyperbolic chains and their monodromy spectrum:\n');
disp(eigenTable);

% Consistency: det = 1, reciprocal eigenvalues, eigenvector residuals.
maxDetError = 0;
maxReciprocalError = 0;
maxTraceError = 0;
maxDirectionResidual = 0;
maxAverageYError = 0;
for i = 1:numel(chains)
    c = chains(i);
    maxAverageYError = max(maxAverageYError, c.averageYError);
    M = c.monodromy;
    maxDetError = max(maxDetError, abs(det(M) - 1));
    maxTraceError = max(maxTraceError, abs(trace(M) - c.trace));
    if isfield(c, 'eigen') && ~isempty(c.eigen)
        e = c.eigen;
        maxReciprocalError = max(maxReciprocalError, ...
            abs(e.stableEigenvalue * e.unstableEigenvalue - 1));
        if isfield(e, 'directionResidual')
            r = e.directionResidual(isfinite(e.directionResidual));
            if ~isempty(r)
                maxDirectionResidual = max(maxDirectionResidual, max(r));
            end
        end
        for j = 1:c.period
            Mj = e.monodromyByPhasePoint(:, :, j);
            maxDetError = max(maxDetError, abs(det(Mj) - 1));
            maxTraceError = max(maxTraceError, abs(trace(Mj) - c.trace));
        end
    end
end

fprintf('\nspectrum consistency over the whole catalogue\n');
fprintf('  chains                                : %d\n', numel(chains));
fprintf('  max |det M - 1| (all phase points)    : %.3e\n', maxDetError);
fprintf('  max |trace M_j - trace M_1|           : %.3e\n', maxTraceError);
fprintf('  max |lambda_s * lambda_u - 1|         : %.3e\n', maxReciprocalError);
fprintf('  max eigenvector residual             : %.3e\n', maxDirectionResidual);
fprintf('  max |mean(y) - m/p| over catalogue    : %.3e\n', maxAverageYError);

try
    fig = so_plot_eigendirections(chains, cfg);
    fprintf('eigendirection figure drawn (handle %g)\n', double(fig.Number));
    fprintf('PART B: PASSED\n');
catch err
    fprintf('PART B: FIGURE FAILED - %s\n', err.message);
    fprintf('%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
end

% ---------------------------------------------------------------- part C
fprintf('\n\n############ PART C: paper-trial timing ############\n');
parts = "census";
if runHotspots
    parts(end + 1) = "hotspots";
end
if runStages
    parts(end + 1) = "stages";
end
try
    diagnostics = so_profile_paper_trial(parts, 1, 1800); %#ok<NASGU>
    fprintf('PART C: COMPLETED\n');
catch err
    fprintf('PART C: FAILED - %s\n', err.message);
    fprintf('%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
end

fprintf('\n===============================================================\n');
fprintf('round 1 diagnostics finished\n');
fprintf('log: %s\n', logFile);
fprintf('===============================================================\n');
diary off;
end
