%RUN_LORENZ_COMPARISON  Table-1 style comparison for the Lorenz system.
%
%   Builds ECR-II and ECR-III models with 1,2,...,Kmax activation regions and
%   reports, for each of them,
%
%       * the off-line training time,
%       * the average reaching time without noise,
%       * the average reaching time with measurement noise,
%
%   together with the uncontrolled system as a reference (the row that plays
%   the role of the "OGY" column of Table 1: no targeting, simply wait until
%   the chaotic trajectory enters the target region by itself).
%
%   The table is printed and written to results/lorenz_comparison.csv.
%
%   Run  startup_ecr  first.  Takes a few minutes.

more off
COORD = 'real';          % 'real' or 'delay'
KLIST = 1:5;

switch COORD
    case 'real',  S = sys_lorenz_real();
    case 'delay', S = sys_lorenz_delay();
end
fprintf('\n=== ECR comparison, %s ===\n%s\n\n', S.name, S.notes);

opt          = ecr_default_options();
opt.nData    = 6000;
opt.trajLen  = 30;
opt.nTrials  = 100;
opt.maxSteps = 150;
opt.verbose  = 0;
if strcmp(COORD,'delay')
    opt.nData   = 3000;
    opt.nTrials = 60;
end
noiseLevel = 0.0008;                 % 0.08 % (Table 1, Lorenz real coord.)

D  = ecr_cache_data(S, setfield(opt, 'verbose', 1));  %#ok<SFLD>
ecr_seed(opt.seed + 100);
H0 = S.init(opt.nTrials);

optN = opt;  optN.noise = noiseLevel;

b0  = ecr_reaching_time(S, [], opt,  H0);
b0n = ecr_reaching_time(S, [], optN, H0);
fprintf('no targeting            : reach %6.2f (noiseless) %6.2f (noisy)\n\n', ...
        b0.meanCens, b0n.meanCens);

rows = {};
fprintf('%-9s %3s %10s %8s %8s %8s %8s\n', ...
        'method','K','train [s]','nets','reach','reach(n)','success');
for K = KLIST
    o = opt;  o.Kmax = K;
    for v = {'ECR-II','ECR-III'}
        variant = v{1};
        m  = ecr_train(S, D, o, variant);
        c  = @(Z) ecr_control(m, Z);
        bc = ecr_reaching_time(S, c, opt,  H0);
        bn = ecr_reaching_time(S, c, optN, H0);
        fprintf('%-9s %3d %10.2f %8d %8.2f %8.2f %8.2f\n', ...
                variant, numel(m.levels)-1, m.info.trainTime, m.info.nNets, ...
                bc.meanCens, bn.meanCens, bc.success);
        rows(end+1,:) = {variant, numel(m.levels)-1, m.info.trainTime, ...
                         m.info.nNets, bc.meanCens, bn.meanCens, bc.success}; %#ok<AGROW>
    end
end

% ---- write csv ----------------------------------------------------------
resdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
if exist(resdir,'dir') ~= 7, mkdir(resdir); end
fn = fullfile(resdir, sprintf('%s_comparison.csv', S.name));
fid = fopen(fn, 'w');
fprintf(fid, 'method,regions,train_time_s,networks,reach_noiseless,reach_noisy,success\n');
fprintf(fid, 'none,0,0,0,%.3f,%.3f,%.3f\n', b0.meanCens, b0n.meanCens, b0.success);
for i = 1:size(rows,1)
    fprintf(fid, '%s,%d,%.3f,%d,%.3f,%.3f,%.3f\n', rows{i,1}, rows{i,2}, ...
            rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7});
end
fclose(fid);
fprintf('\nwritten %s\n', fn);
