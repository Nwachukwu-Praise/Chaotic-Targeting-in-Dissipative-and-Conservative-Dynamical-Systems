function result = main_schroer_ott_targeting()
%MAIN_SCHROER_OTT_TARGETING First-light Schroer-Ott standard-map run.
%
% This is intentionally a script-like function with a single central
% configuration.  It does not create an .mlx file.
cfg = schroer_ott_default_config();
rng(cfg.randomSeed);
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end

fprintf('Schroer-Ott first-light targeting on S^1 x R\n');
fprintf('k: %.12g\n', cfg.k);
fprintf('control amplitude: %.12g\n', cfg.controlAmplitude);
fprintf('maxPeriod: %d\n', cfg.orbit.maxPeriod);

catalogue = so_enumerate_periodic_orbits(cfg);
disp('Orbit catalogue sorted by omega:');
disp(catalogue.byOmega);
disp('Diagnostic Figure-3 fractions:');
disp(catalogue.diagnosticFractions);

route = so_construct_first_light_route(catalogue, cfg);
fprintf('route rotation numbers: %s\n', mat2str(route.rotationNumbers, 6));
if ~isempty(route.ambiguities) && height(route.ambiguities) > 0
    disp('Route ambiguities retained for diagnostics:');
    disp(route.ambiguities);
end

result = so_multistage_targeting(cfg, catalogue, route);
result.diagnosticManifolds = so_compute_diagnostic_manifolds(route, cfg);
so_print_result_summary(result);
so_save_outputs(result);
so_plot_results(result);
end
