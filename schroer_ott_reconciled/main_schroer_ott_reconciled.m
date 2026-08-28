function study = main_schroer_ott_reconciled(varargin)
%MAIN_SCHROER_OTT_RECONCILED Full reconciled Schroer-Ott demonstration.
%
%   study = MAIN_SCHROER_OTT_RECONCILED()
%   study = MAIN_SCHROER_OTT_RECONCILED('Ensemble', false)
%   study = MAIN_SCHROER_OTT_RECONCILED('EnsembleTrials', 5, 'Realisations', 20)
%
% Runs, in order:
%   1. the three source/target case studies (diagonal, horizontal, vertical)
%   2. the 50-source ensemble on the diagonal geometry
%   3. the additive-noise sweep on the diagonal solution
%
% Everything is written under ./outputs.  Rough timings, from the 23 s the
% diagonal case takes as a yardstick:
%   diagonal     under a minute
%   vertical     a few minutes  (three stages, more switch probes)
%   horizontal   the slow one   (no proxy shortens the final leg, so the
%                               forward-backward search runs out to tau = 27
%                               against curves of thousands of points)
%   ensemble     20-60 minutes  (checkpointed; safe to interrupt and resume)
%   noise        under a minute
%
% For a quick smoke test that touches every code path:
%   main_schroer_ott_reconciled('EnsembleTrials', 3, 'Realisations', 10)

opts = so_parse_options(varargin, struct( ...
    'Cases', {{'diagonal', 'horizontal', 'vertical'}}, ...
    'Ensemble', true, 'EnsembleCase', 'diagonal', 'EnsembleTrials', [], ...
    'Uncontrolled', true, 'Noise', true, 'NoiseCase', 'diagonal', ...
    'Realisations', [], 'Plot', true, 'Save', true, 'Verbose', true));

totalTimer = tic;
fprintf('Schroer-Ott pass targeting, standard map on S^1 x R\n');
fprintf('reconciled build: fundamental targeting engine + phase-portrait\n');
fprintf('figures + three-geometry, ensemble and noise validations\n');

study.cases = so_run_all_cases('Cases', opts.Cases, 'Plot', opts.Plot, ...
    'Save', opts.Save, 'Verbose', opts.Verbose);

study.ensemble = [];
if opts.Ensemble
    args = {'Case', opts.EnsembleCase, 'Uncontrolled', opts.Uncontrolled, ...
        'Plot', opts.Plot, 'Save', opts.Save, 'Verbose', opts.Verbose};
    if ~isempty(opts.EnsembleTrials)
        args = [args, {'TrialCount', opts.EnsembleTrials}];
    end
    study.ensemble = so_run_source_ensemble(args{:});
end

study.noise = [];
if opts.Noise
    idx = find(strcmpi(opts.Cases, opts.NoiseCase), 1);
    if isempty(idx)
        warning('SchroerOtt:NoiseCaseNotRun', ...
            'Noise case "%s" was not among the cases run; skipping noise validation.', ...
            opts.NoiseCase);
    else
        args = {'Plot', opts.Plot, 'Save', opts.Save, 'Verbose', opts.Verbose};
        if ~isempty(opts.Realisations)
            args = [args, {'Realisations', opts.Realisations}];
        end
        study.noise = so_run_noise_validation(study.cases.results{idx}, args{:});
    end
end

study.runtimeSeconds = toc(totalTimer);
fprintf('\nComplete in %.1f s.  Outputs under %s\n', study.runtimeSeconds, ...
    fullfile(pwd, 'outputs'));

if opts.Save
    outDir = fullfile(pwd, 'outputs');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    save(fullfile(outDir, 'reconciled_study.mat'), 'study', '-v7.3');
end
end
