%DEMO_OPTIONS_STUDY  Sensitivity of ECR-III to the details the paper leaves open.
%
%   Three things are not fixed by Iplikci & Denizhan (2003):
%     * how the "Normalized Covariance Matrix" of Eq. (8) is normalised
%       (OPT.covNorm),
%     * whether a cluster may be rejected when the current state is far
%       outside the cloud it was trained on (OPT.nmdGate),
%     * which cluster wins when several regions could claim the state
%       (OPT.selection: nearest cluster overall, or lowest region index).
%
%   Every configuration is trained NREP times (the network centres are placed
%   by a randomised k-means, so a single run is a noisy measurement) and the
%   reaching time is averaged over the repetitions, always from the same set
%   of initial conditions.  Results go to results/lorenz_real_options.csv and
%   justify the defaults in ECR_DEFAULT_OPTIONS.
%
%   Run  startup_ecr  first.  Takes ~20 minutes.

more off
NREP = 3;

S = sys_lorenz_real();

opt          = ecr_default_options();
opt.nData    = 6000;
opt.trajLen  = 30;
opt.nTrials  = 100;
opt.maxSteps = 150;
opt.Kmax     = 5;
opt.verbose  = 0;

D = ecr_cache_data(S, setfield(opt,'verbose',1));   %#ok<SFLD>

ecr_seed(12345);                    % identical initial conditions everywhere
H0 = S.init(opt.nTrials);

b0 = ecr_reaching_time(S, [], opt, H0);
fprintf('\nno targeting: %.2f steps (success %.2f)\n\n', b0.meanCens, b0.success);

norms = {'det','plain','diag','trace'};
gates = [1.5 2.0 Inf];

rows = {};
fprintf('--- part A: normalisation of Eq. (8) and the NMD gate (selection = nmd)\n');
fprintf('%-6s %6s | %14s %8s %8s %6s\n', 'covN','gate','reach +- std','success','hold','nets');
best = struct('reach', Inf);
for i = 1:numel(norms)
    for g = gates
        o = opt;  o.covNorm = norms{i};  o.nmdGate = g;  o.selection = 'nmd';
        [r, sc, hd, nn] = ecr_repeat_config(S, D, o, H0, NREP);
        fprintf('%-6s %6.1f | %7.2f +-%5.2f %8.2f %8.2f %6.1f\n', ...
                norms{i}, g, mean(r), std(r), mean(sc), mean(hd), mean(nn));
        rows(end+1,:) = {norms{i}, g, 'nmd', mean(r), std(r), mean(sc), mean(hd), mean(nn)}; %#ok<AGROW>
        if mean(r) < best.reach
            best = struct('reach', mean(r), 'covNorm', norms{i}, 'gate', g);
        end
    end
end

fprintf('\n--- part B: region selection rule, at covNorm = %s, gate = %.1f\n', ...
        best.covNorm, best.gate);
for sel = {'nmd','level'}
    o = opt;  o.covNorm = best.covNorm;  o.nmdGate = best.gate;  o.selection = sel{1};
    [r, sc, hd, nn] = ecr_repeat_config(S, D, o, H0, NREP);
    fprintf('%-6s        | %7.2f +-%5.2f %8.2f %8.2f %6.1f\n', ...
            sel{1}, mean(r), std(r), mean(sc), mean(hd), mean(nn));
    rows(end+1,:) = {best.covNorm, best.gate, sel{1}, mean(r), std(r), ...
                     mean(sc), mean(hd), mean(nn)}; %#ok<AGROW>
end

resdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
if exist(resdir,'dir') ~= 7, mkdir(resdir); end
fn  = fullfile(resdir, 'lorenz_real_options.csv');
fid = fopen(fn,'w');
fprintf(fid, 'covNorm,nmdGate,selection,reach_mean,reach_std,success,hold,networks\n');
fprintf(fid, 'none,0,none,%.3f,0,%.3f,0,0\n', b0.meanCens, b0.success);
for i = 1:size(rows,1)
    fprintf(fid, '%s,%g,%s,%.3f,%.3f,%.3f,%.3f,%.1f\n', rows{i,1}, rows{i,2}, ...
            rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7}, rows{i,8});
end
fclose(fid);
fprintf('\nwritten %s\n', fn);
