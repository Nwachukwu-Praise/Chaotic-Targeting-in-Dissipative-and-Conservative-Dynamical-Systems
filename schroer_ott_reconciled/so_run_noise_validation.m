function noise = so_run_noise_validation(result, varargin)
%SO_RUN_NOISE_VALIDATION Robustness of a solved case to additive noise.
%
%   noise = SO_RUN_NOISE_VALIDATION(result)
%   noise = SO_RUN_NOISE_VALIDATION(result, 'Realisations', 200)
%
% Schroer and Ott: "We studied targeting the standard map under the
% influence of small additive noise with amplitudes ranging up to 1e-2 x
% delta at each iteration.  Without any change to our targeting scheme, all
% trajectories reached the target."  Here delta is the full length of the
% admissible control segment, 2 * cfg.controlAmplitude = 0.006, so the
% headline noise level is sigma = 6e-5 per iteration.
%
% The test sweeps several multiples of delta rather than only the published
% one, because a single level that happens to pass proves less than a curve
% showing where the scheme starts to break.  sigma = 0 is included as a
% control: it must reproduce the deterministic endpoint exactly, which is
% also a check that the replay itself is faithful.
%
% Name/value options:
%   'SigmaFactors'  multiples of delta to sweep   (default from cfg.noise)
%   'Realisations'  independent draws per level   (default from cfg.noise)
%   'Seed'          RNG seed                      (default from cfg.noise)
%   'Save'          write MAT/CSV                 (default true)
%   'Plot'          draw the figure               (default true)
%   'Verbose'       print the summary             (default true)

cfg = result.configuration;
opts = so_parse_options(varargin, struct( ...
    'SigmaFactors', cfg.noise.sigmaFactors, ...
    'Realisations', cfg.noise.realisationsPerLevel, ...
    'Seed', cfg.noise.randomSeed, ...
    'Save', true, 'Plot', true, 'Verbose', true));

if ~result.targetReached
    error('SchroerOtt:NoiseOnFailedRun', ...
        'The deterministic run did not reach the target; there is no schedule to perturb.');
end

delta = 2 * cfg.controlAmplitude;
factors = opts.SigmaFactors(:).';
totalIterations = sum([result.executionSegments.iterations]);

rng(opts.Seed);
% One noise bank shared across levels: the same realisation is then scaled
% to each sigma, so differences between levels are the amplitude alone and
% not a different random draw.
bank = randn(opts.Realisations, max(1, totalIterations));

rows = {};
records = struct('factor', {}, 'sigma', {}, 'replays', {});
for f = 1:numel(factors)
    sigma = factors(f) * delta;
    replays = cell(1, opts.Realisations);
    contained = false(1, opts.Realisations);
    displacement = zeros(1, opts.Realisations);
    for r = 1:opts.Realisations
        rep = so_replay_with_noise(result, sigma, bank(r, :));
        replays{r} = rep;
        contained(r) = rep.targetContained;
        displacement(r) = rep.endpointDisplacement;
    end
    rows(end + 1, :) = {factors(f), sigma, opts.Realisations, sum(contained), ...
        mean(contained), min(displacement), median(displacement), ...
        max(displacement), std(displacement)}; %#ok<AGROW>
    records(end + 1) = struct('factor', factors(f), 'sigma', sigma, ...
        'replays', {replays}); %#ok<AGROW>
    if opts.Verbose
        fprintf('  sigma = %8.3g (%6.4g x delta): %3d/%3d contained, median endpoint shift %.3g\n', ...
            sigma, factors(f), sum(contained), opts.Realisations, median(displacement));
    end
end

levels = cell2table(rows, 'VariableNames', {'sigmaOverDelta','sigma', ...
    'realisations','containedCount','containedFraction','minDisplacement', ...
    'medianDisplacement','maxDisplacement','stdDisplacement'});

zeroRow = find(factors == 0, 1);
if isempty(zeroRow)
    zeroNoiseExact = NaN;
else
    zeroNoiseExact = levels.maxDisplacement(zeroRow);
end

noise.caseName = cfg.caseName;
noise.delta = delta;
noise.headlineFactor = cfg.noise.headlineFactor;
noise.totalIterations = totalIterations;
noise.randomSeed = opts.Seed;
noise.noiseApplication = "eta_n added to y_{n+1} before x_{n+1} = x_n + y_{n+1}";
noise.schedule = "open loop: deterministic controls replayed unchanged";
noise.levels = levels;
noise.records = records;
noise.zeroNoiseMaxDisplacement = zeroNoiseExact;
noise.deterministicFinalState = result.finalState;

if opts.Verbose
    fprintf('\n================ NOISE VALIDATION (%s) ================\n', char(cfg.caseName));
    fprintf('delta = %g, schedule replayed open loop over %d iterations\n', delta, totalIterations);
    if isfinite(zeroNoiseExact)
        fprintf('zero-noise replay reproduces the deterministic endpoint to %.3g\n', zeroNoiseExact);
    end
    disp(levels);
end

if opts.Save
    outDir = fullfile(cfg.outputDirectory, 'noise');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    save(fullfile(outDir, 'noise_validation.mat'), 'noise', '-v7.3');
    writetable(levels, fullfile(outDir, 'noise_levels.csv'));
end
if opts.Plot
    so_plot_noise(result, noise, opts.Save);
end
end
