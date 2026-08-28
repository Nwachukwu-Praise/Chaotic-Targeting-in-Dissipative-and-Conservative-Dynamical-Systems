function [cfg, paper] = so_paper_comparison_config()
%SO_PAPER_COMPARISON_CONFIG Estimated Figure-3 Schroer-Ott geometry.
%
% The paper does not tabulate rectangle coordinates.  These values are
% digitized by reading the rendered Figure 3 axes and gray rectangles.
cfg = schroer_ott_default_config();
cfg.sourceRectangle = struct('xMin', 0.455, 'xMax', 0.555, ...
    'yMin', 0.000, 'yMax', 0.045, 'id', 'paper-source');
cfg.targetRectangle = struct('xMin', 0.455, 'xMax', 0.555, ...
    'yMin', 0.945, 'yMax', 1.000, 'id', 'final-target');
cfg.maxForwardIterations = 10;
cfg.maxBackwardIterations = 8;
cfg.maxTotalTransferTime = 18;
cfg.curve.initialControlSamples = 17;
cfg.curve.initialBoundarySamples = 64;
cfg.outputDirectory = fullfile(pwd, 'outputs', 'paper_comparison');
cfg.figureDirectory = fullfile(cfg.outputDirectory, 'figures');

paper.source = cfg.sourceRectangle;
paper.target = cfg.targetRectangle;
paper.digitisationMethod = "visual reading of rendered Schroer-Ott Figure 3 page against axes";
paper.coordinateUncertainty = 0.015;
paper.expectedRotationNumbers = [1/4 1/3 2/5 1/2 3/5 2/3 3/4];
paper.publishedControlledRange = [125 132];
paper.publishedExampleIterations = 125;
paper.maxControlAmplitude = cfg.controlAmplitude;
end
