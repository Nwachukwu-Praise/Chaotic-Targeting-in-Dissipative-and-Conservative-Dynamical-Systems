function paperCase = so_run_paper_comparison_case(saveOutputs)
%SO_RUN_PAPER_COMPARISON_CASE Attempt the k=1.25 Figure-3 comparison case.
if nargin < 1
    saveOutputs = true;
end
[cfg, paper] = so_paper_comparison_config();
cfg.saveFigures = saveOutputs;
if saveOutputs
    ensure_output_dirs(cfg);
end
rng(cfg.randomSeed);
timer = tic;
catalogue = so_enumerate_periodic_orbits(cfg);
route = so_construct_omega_route(catalogue, cfg, paper.expectedRotationNumbers, ...
    'estimated Figure-3 seven-resonance route');
try
    result = so_multistage_targeting(cfg, catalogue, route);
    if result.targetReached
        classification = classify_paper_reproduction(result, paper);
    else
        classification = "failed numerical reproduction";
    end
catch err
    result = [];
    classification = "failed numerical reproduction";
    paperCase.errorIdentifier = string(err.identifier);
    paperCase.errorMessage = string(err.message);
end
paperCase.configuration = cfg;
paperCase.paperMetadata = paper;
paperCase.catalogue = catalogue;
paperCase.route = route;
paperCase.result = result;
paperCase.classification = classification;
paperCase.runtimeSeconds = toc(timer);
if saveOutputs
    save(fullfile(cfg.outputDirectory, 'paper_comparison_case.mat'), 'paperCase', '-v7.3');
end
end

function classification = classify_paper_reproduction(result, paper)
if result.targetReached && result.totalExecutedIterations >= paper.publishedControlledRange(1) && ...
        result.totalExecutedIterations <= paper.publishedControlledRange(2) && ...
        numel(result.route.rotationNumbers) == numel(paper.expectedRotationNumbers)
    classification = "partial reproduction";
elseif result.targetReached
    classification = "qualitative reproduction";
else
    classification = "failed numerical reproduction";
end
end

function ensure_output_dirs(cfg)
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
if ~exist(cfg.figureDirectory, 'dir')
    mkdir(cfg.figureDirectory);
end
end
