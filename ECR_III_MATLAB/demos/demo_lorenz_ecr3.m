%DEMO_LORENZ_ECR3  ECR-III targeting of the Lorenz system (real coordinates).
%
%   Reproduces the "Lorenz system (real coord.)" experiment of Iplikci &
%   Denizhan (2003):  the flow is reduced to the Poincare map on the surface
%   y = 8.4853, the target is the period-1 unstable orbit at
%   z* = [14.2387 39.7934] (Table 2) and the three control parameters
%   sigma, rho, beta may be perturbed by at most [0.30 0.84 0.08].
%
%   The script
%     1. collects (z_n, p_n, z_{n+1}) triples with random parameters in Pi,
%     2. builds the extended control regions T_0..T_K and their clusters,
%     3. trains one RBF network per cluster,
%     4. measures the average reaching time and compares it with the
%        uncontrolled system and with ECR-II,
%     5. saves figures into results/.
%
%   Run  startup_ecr  first.

more off
fprintf('\n=== ECR-III on the Lorenz system (real coordinates) ===\n\n');

%% ---- 1. system and options --------------------------------------------
S = sys_lorenz_real();
fprintf('%s\n\n', S.notes);

opt          = ecr_default_options();
opt.nData    = 6000;      % raw (z,p,z') triples
opt.trajLen  = 30;
opt.Kmax     = 5;         % T_1 ... T_5   (Table 1 uses 1..6 regions)
opt.nTrials  = 100;       % initial conditions used for the reaching time
opt.maxSteps = 150;
opt.seed     = 1;

%% ---- 2. data -----------------------------------------------------------
D = ecr_cache_data(S, opt);

%% ---- 3. training -------------------------------------------------------
m3 = ecr_train(S, D, opt, 'ECR-III');
m2 = ecr_train(S, D, opt, 'ECR-II');

fprintf('\n');
ecr_summary(m3);
fprintf('\n');

%% ---- 4. reaching times -------------------------------------------------
ecr_seed(opt.seed + 100);                 % reproducible test set
H0 = S.init(opt.nTrials);                 % same initial conditions for all

optN       = opt;                         % noisy variant of Table 1
optN.noise = 0.0008;                      % 0.08 % of the measured RMS

fprintf('\n--- average reaching time over %d initial conditions ---\n', opt.nTrials);
res = struct();
for noisy = [false true]
    o = opt; if noisy, o = optN; end
    b0 = ecr_reaching_time(S, [],                        o, H0);
    b2 = ecr_reaching_time(S, @(Z) ecr_control(m2,Z),    o, H0);
    b3 = ecr_reaching_time(S, @(Z) ecr_control(m3,Z),    o, H0);
    tag = 'noiseless'; if noisy, tag = '0.08% noise'; end
    fprintf('  %-12s  no control %6.2f | ECR-II %6.2f | ECR-III %6.2f\n', ...
            tag, b0.meanCens, b2.meanCens, b3.meanCens);
    fprintf('  %-12s  success     %6.2f | %6.2f | %6.2f   (hold %.2f / %.2f)\n', ...
            '', b0.success, b2.success, b3.success, b2.hold, b3.hold);
    if noisy, res.noisy = {b0,b2,b3}; else, res.clean = {b0,b2,b3}; end
end

%% ---- 5. one closed-loop run and the figures ---------------------------
resdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
if exist(resdir,'dir') ~= 7, mkdir(resdir); end

h0 = S.init(1);
R3 = ecr_simulate(S, @(z) ecr_control(m3,z), h0, opt);
R0 = ecr_simulate(S, [],                     h0, opt);
fprintf('\n  example run: ECR-III reached S0 in %s steps, uncontrolled %s steps\n', ...
        num2str(R3.reach), num2str(R0.reach));

ecr_plot_clusters(m3, D, [1 2], fullfile(resdir,'lorenz_real_clusters.png'));
ecr_plot_run(S, R3,              fullfile(resdir,'lorenz_real_run.png'));
fprintf('  figures written to %s\n', resdir);

ecr_save_model(fullfile(resdir,'lorenz_real_models.mat'), m2, m3);
fprintf('\n=== done ===\n');
