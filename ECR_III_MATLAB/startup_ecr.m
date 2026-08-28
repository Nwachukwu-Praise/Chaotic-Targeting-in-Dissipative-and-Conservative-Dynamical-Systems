function startup_ecr()
%STARTUP_ECR  Put the ECR toolbox on the MATLAB/Octave path.
%
%   Run this once per session from the ECR_III_MATLAB folder:
%
%       >> startup_ecr
%       >> demo_lorenz_ecr3          % main demonstration
%
%   Works in MATLAB and in GNU Octave (no toolboxes are required).

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'core'));
addpath(fullfile(here, 'systems'));
addpath(fullfile(here, 'demos'));
addpath(fullfile(here, 'tools'));
addpath(fullfile(here, 'tests'));
addpath(fullfile(here, 'tests', 'fixtures'));
fprintf('ECR-III toolbox on path (%s)\n', here);
end
