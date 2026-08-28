function [cfg, meta] = so_case_config(caseName)
%SO_CASE_CONFIG Geometry for one of the three source/target case studies.
%
%   [cfg, meta] = SO_CASE_CONFIG('diagonal')     % the original demonstration
%   [cfg, meta] = SO_CASE_CONFIG('horizontal')   % same transport band
%   [cfg, meta] = SO_CASE_CONFIG('vertical')     % same x-band, y-transport
%
% The three geometries answer the question the single diagonal run cannot:
% does the pass-targeting method work because the method works, or because
% that one source/target pair happened to be favourable?
%
% All three rectangles are 0.04-0.05 on a side and all six centres lie in the
% connected chaotic region of the k = 1.25 standard map (finite-time Lyapunov
% exponent 0.18-0.22 per iterate; a regular island would give ~0).  Being a
% mixed phase space, every rectangle of this size also contains small regular
% fragments -- so does the original diagonal target -- which is why
% containment is tested on the trajectory, never assumed from the rectangle.
%
% The iteration budgets differ per case and are not arbitrary.  Equation (1)
% of the paper gives tau ~ lambda1^-1 ln(L/delta) + |lambda2|^-1 ln(L/eps_t);
% with lambda ~ 0.2 per iterate, delta = 0.006 and eps_t = 0.04 that is
% roughly 40 iterates for a connection made in one step, and far less when a
% resonance proxy sits near the target.  The budgets below were set from a
% forward-backward reachability sweep over each geometry, and are the
% smallest round values that leave headroom above the observed minimum.

if nargin < 1 || isempty(caseName)
    caseName = 'diagonal';
end
caseName = lower(char(caseName));

cfg = so_reconciled_config();

switch caseName
    case 'diagonal'
        cfg.caseName = "diagonal";
        cfg.caseLabel = "Diagonal: source low-right, target high-left";
        cfg.sourceRectangle = struct('xMin', 0.64, 'xMax', 0.69, ...
            'yMin', 0.48, 'yMax', 0.52, 'id', 'source');
        cfg.targetRectangle = struct('xMin', 0.21, 'xMax', 0.25, ...
            'yMin', 0.71, 'yMax', 0.75, 'id', 'final-target');
        cfg.maxForwardIterations = 11;
        cfg.maxBackwardIterations = 8;
        cfg.maxTotalTransferTime = 19;
        meta.expectedRotationNumbers = [3/5, 2/3];
        meta.note = ['Reference case, identical to the fundamental-folder ' ...
            'demonstration.  Reproduces its result exactly: two proxy ' ...
            'resonances, the first skipped, 9 executed iterations.'];

    case 'horizontal'
        cfg.caseName = "horizontal";
        cfg.caseLabel = "Horizontal: same transport band, half a period apart in x";
        cfg.sourceRectangle = struct('xMin', 0.315, 'xMax', 0.355, ...
            'yMin', 0.28, 'yMax', 0.32, 'id', 'source');
        cfg.targetRectangle = struct('xMin', 0.815, 'xMax', 0.855, ...
            'yMin', 0.28, 'yMax', 0.32, 'id', 'final-target');
        % Measured: this case resolves at tau = 27, split (nForward, nBackward)
        % = (11, 16), control -0.000387497830700472, in about 15 minutes.  It
        % is the slowest of the three because nothing shortens the final leg:
        % with no proxy near the target the forward-backward search has to run
        % all the way out to tau = 27.
        %
        % Do not trim these budgets to the winning split.  The forward images
        % beyond n = 13 are large (9 330 points at n = 13, then 15 689, 25 935
        % and 42 798) and are returned unresolved, which is tempting to cut --
        % but the accepted connection uses nBackward = 16, so capping the
        % backward budget lower forces a different, larger-control solution at
        % the same tau.  Same transfer time, worse control, for a speed-up
        % that is not worth changing the answer.
        cfg.maxForwardIterations = 16;
        cfg.maxBackwardIterations = 16;
        cfg.maxTotalTransferTime = 30;
        meta.expectedRotationNumbers = zeros(1, 0);
        meta.note = ['Source and target share the y-band, so no resonance ' ...
            'separates them and the rotation bracket is empty.  The route ' ...
            'has no intermediate targets and the method reduces to a single ' ...
            'Shinbrot forward-backward step.  This is the degenerate limit ' ...
            'of pass targeting and is the case that most directly tests ' ...
            'the forward-backward core on its own.'];

    case 'vertical'
        cfg.caseName = "vertical";
        cfg.caseLabel = "Vertical: same x-band, transport across three resonances";
        cfg.sourceRectangle = struct('xMin', 0.425, 'xMax', 0.465, ...
            'yMin', 0.40, 'yMax', 0.44, 'id', 'source');
        cfg.targetRectangle = struct('xMin', 0.425, 'xMax', 0.465, ...
            'yMin', 0.70, 'yMax', 0.74, 'id', 'final-target');
        cfg.maxForwardIterations = 16;
        cfg.maxBackwardIterations = 14;
        cfg.maxTotalTransferTime = 30;
        meta.expectedRotationNumbers = [1/2, 3/5, 2/3];
        meta.note = ['Pure y-transport at fixed x.  The bracket now holds ' ...
            'three direct-hyperbolic resonances instead of two, so this is ' ...
            'the case that exercises multistage switching hardest, and the ' ...
            'one closest in spirit to Figure 3 of the paper.'];

    otherwise
        error('SchroerOtt:UnknownCase', ...
            'Unknown case "%s".  Use diagonal, horizontal or vertical.', caseName);
end

meta.caseName = cfg.caseName;
meta.caseLabel = cfg.caseLabel;
meta.sourceRectangle = cfg.sourceRectangle;
meta.targetRectangle = cfg.targetRectangle;
meta.sourceCentre = so_rectangle_center(cfg.sourceRectangle);
meta.targetCentre = so_rectangle_center(cfg.targetRectangle);
meta.iterationBudget = [cfg.maxForwardIterations, cfg.maxBackwardIterations, ...
    cfg.maxTotalTransferTime];

cfg.outputDirectory = fullfile(pwd, 'outputs', char(cfg.caseName));
cfg.figureDirectory = fullfile(cfg.outputDirectory, 'figures');
end
