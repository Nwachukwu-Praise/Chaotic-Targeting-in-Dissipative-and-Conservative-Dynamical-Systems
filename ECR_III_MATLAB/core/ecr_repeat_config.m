function [r, sc, hd, nn] = ecr_repeat_config(S, D, opt, H0, nrep, variant)
%ECR_REPEAT_CONFIG  Train and benchmark one configuration several times.
%
%   [R,SC,HD,NN] = ECR_REPEAT_CONFIG(S,D,OPT,H0,NREP,VARIANT)
%
%   The RBF centres are placed by a randomised k-means, so a single training
%   run is a noisy measurement of a configuration's quality.  This helper
%   repeats training NREP times (different network seeds, identical data and
%   identical initial conditions H0) and returns the vectors of
%
%     R   average reaching time         SC  capture rate
%     HD  retention after capture       NN  number of networks
%
%   Used by DEMO_OPTIONS_STUDY.

if nargin < 5 || isempty(nrep),    nrep = 3; end
if nargin < 6 || isempty(variant), variant = 'ECR-III'; end

r = zeros(1,nrep); sc = r; hd = r; nn = r;
for k = 1:nrep
    ecr_seed(1000 + k);
    m = ecr_train(S, D, opt, variant);
    b = ecr_reaching_time(S, @(Z) ecr_control(m,Z), opt, H0);
    r(k)  = b.meanCens;
    sc(k) = b.success;
    hd(k) = b.hold;
    nn(k) = m.info.nNets;
end
end
