function [result, meta] = so_run_case(caseName, varargin)
%SO_RUN_CASE Run one source/target case study end to end.
%
%   result = SO_RUN_CASE('diagonal')
%   result = SO_RUN_CASE('vertical', 'Catalogue', cat, 'Plot', false)
%
% Name/value options:
%   'Catalogue'  pre-computed orbit catalogue (avoids re-enumerating)
%   'Plot'       produce figures                        (default true)
%   'Save'       write MAT/CSV outputs                  (default true)
%   'Verbose'    print the console summary              (default true)
%   'Manifolds'  compute a posteriori manifold branches (default true)

opts = so_parse_options(varargin, struct( ...
    'Catalogue', [], 'Plot', true, 'Save', true, 'Verbose', true, 'Manifolds', true));

[cfg, meta] = so_case_config(caseName);
cfg.saveFigures = opts.Plot && opts.Save;
cfg.verbose = opts.Verbose;
rng(cfg.randomSeed);

if opts.Save
    if ~exist(cfg.outputDirectory, 'dir'), mkdir(cfg.outputDirectory); end
    if ~exist(cfg.figureDirectory, 'dir'), mkdir(cfg.figureDirectory); end
end

if opts.Verbose
    fprintf('\n================================================================\n');
    fprintf('CASE: %s\n', char(cfg.caseLabel));
    fprintf('  source rectangle  x=[%.4g, %.4g]  y=[%.4g, %.4g]\n', ...
        cfg.sourceRectangle.xMin, cfg.sourceRectangle.xMax, ...
        cfg.sourceRectangle.yMin, cfg.sourceRectangle.yMax);
    fprintf('  target rectangle  x=[%.4g, %.4g]  y=[%.4g, %.4g]\n', ...
        cfg.targetRectangle.xMin, cfg.targetRectangle.xMax, ...
        cfg.targetRectangle.yMin, cfg.targetRectangle.yMax);
    fprintf('  budgets  nForward<=%d  nBackward<=%d  tau<=%d\n', ...
        cfg.maxForwardIterations, cfg.maxBackwardIterations, cfg.maxTotalTransferTime);
    fprintf('================================================================\n');
end

if isempty(opts.Catalogue)
    catalogue = so_enumerate_periodic_orbits(cfg);
else
    catalogue = opts.Catalogue;
end

route = so_construct_route(catalogue, cfg);
if opts.Verbose
    if route.bracketEmpty
        fprintf('route: empty rotation bracket over y in (%.4g, %.4g)\n', ...
            route.bracket(1), route.bracket(2));
        fprintf('       no intermediate resonance -- single forward-backward step\n');
    else
        fprintf('route rotation numbers: %s\n', mat2str(route.rotationNumbers, 6));
    end
end

result = so_multistage_targeting(cfg, catalogue, route);
if opts.Manifolds && ~route.bracketEmpty
    result.diagnosticManifolds = so_compute_diagnostic_manifolds(route, cfg);
end
result.caseMetadata = meta;

if opts.Verbose
    so_print_result_summary(result);
end
if opts.Save
    so_save_outputs(result);
end
if opts.Plot
    so_plot_results(result);
end
end
