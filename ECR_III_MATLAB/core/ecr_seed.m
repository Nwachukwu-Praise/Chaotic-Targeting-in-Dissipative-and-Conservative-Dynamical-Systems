function ecr_seed(s)
%ECR_SEED  Seed the random number generator (MATLAB and Octave compatible).
if exist('rng', 'builtin') == 5 || exist('rng', 'file') == 2
    rng(s, 'twister');
else
    rand('state', s);   %#ok<RAND>
    randn('state', s);  %#ok<RAND>
end
end
