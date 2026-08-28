function study = so_run_all_cases(varargin)
%SO_RUN_ALL_CASES Run the diagonal, horizontal and vertical case studies.
%
%   study = SO_RUN_ALL_CASES()
%   study = SO_RUN_ALL_CASES('Cases', {'diagonal','vertical'}, 'Plot', false)
%
% The point of running three geometries rather than one is falsifiability.
% A single diagonal transfer that succeeds tells you very little: the source
% and target might simply have been well placed.  Three transfers whose
% routes differ in kind -- two resonances, none at all, and three -- test the
% same engine against qualitatively different transport problems:
%
%   diagonal    two resonances crossed, the first one skipped by the
%               switching rule; the reference result
%   horizontal  no resonance between source and target, so the rotation
%               bracket is empty and pass targeting must degenerate
%               gracefully to a plain forward-backward step
%   vertical    pure y-transport across three resonances, the hardest
%               multistage case and the closest to the paper's Figure 3
%
% The orbit catalogue is enumerated once and shared, so the three cases are
% guaranteed to be working from identical periodic-orbit data.

opts = so_parse_options(varargin, struct( ...
    'Cases', {{'diagonal', 'horizontal', 'vertical'}}, ...
    'Plot', true, 'Save', true, 'Verbose', true, 'Manifolds', true));

cases = opts.Cases;
if ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end

timer = tic;
baseCfg = so_reconciled_config();
if opts.Verbose
    fprintf('Enumerating periodic orbits once for all cases...\n');
end
catalogue = so_enumerate_periodic_orbits(baseCfg);
if opts.Verbose
    fprintf('  %d chains, %d direct-hyperbolic\n', numel(catalogue.chains), ...
        sum(strcmp({catalogue.chains.classification}, 'direct-hyperbolic')));
end

results = cell(1, numel(cases));
metas = cell(1, numel(cases));
for i = 1:numel(cases)
    [results{i}, metas{i}] = so_run_case(cases{i}, 'Catalogue', catalogue, ...
        'Plot', opts.Plot, 'Save', opts.Save, 'Verbose', opts.Verbose, ...
        'Manifolds', opts.Manifolds);
end

comparison = so_case_comparison_table(results);

study.cases = cases;
study.catalogue = catalogue;
study.results = results;
study.metadata = metas;
study.comparison = comparison;
study.runtimeSeconds = toc(timer);

if opts.Verbose
    fprintf('\n================ CASE COMPARISON ================\n');
    disp(comparison);
end

outDir = fullfile(pwd, 'outputs');
if opts.Save
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    writetable(comparison, fullfile(outDir, 'case_comparison.csv'));
    save(fullfile(outDir, 'case_study.mat'), 'study', '-v7.3');
end
if opts.Plot
    so_plot_case_comparison(results, opts.Save);
end
end
