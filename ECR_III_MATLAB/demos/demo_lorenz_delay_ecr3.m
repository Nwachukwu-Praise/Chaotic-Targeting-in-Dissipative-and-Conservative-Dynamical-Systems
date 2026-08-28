%DEMO_LORENZ_DELAY_ECR3  ECR-III on the Lorenz system in delay coordinates.
%
%   Same experiment as DEMO_LORENZ_ECR3 but the controller only sees a scalar
%   measurement of the flow, embedded with the delay T = 100 ms into
%
%       z_n = [s(t_n) s(t_n-T) s(t_n-2T)]      (Table 2, "delay coord.")
%
%   This is the practically relevant case: no state observer, no model, only
%   a measured signal sampled at the crossings of the surface of section.
%
%   The run is smaller than the real-coordinate demo because every map step
%   needs the flow to be integrated up to the next crossing while a delay
%   history is recorded.
%
%   Run  startup_ecr  first.

more off
fprintf('\n=== ECR-III on the Lorenz system (delay coordinates) ===\n\n');

S = sys_lorenz_delay();
fprintf('%s\n\n', S.notes);

opt          = ecr_default_options();
opt.nData    = 3000;
opt.trajLen  = 25;
opt.Kmax     = 4;
opt.nTrials  = 60;
opt.maxSteps = 150;
opt.seed     = 1;

D  = ecr_cache_data(S, opt);
m3 = ecr_train(S, D, opt, 'ECR-III');
m2 = ecr_train(S, D, opt, 'ECR-II');
fprintf('\n'); ecr_summary(m3); fprintf('\n');

ecr_seed(opt.seed + 100);
H0 = S.init(opt.nTrials);
b0 = ecr_reaching_time(S, [],                     opt, H0);
b2 = ecr_reaching_time(S, @(Z) ecr_control(m2,Z), opt, H0);
b3 = ecr_reaching_time(S, @(Z) ecr_control(m3,Z), opt, H0);

fprintf('\n--- average reaching time over %d initial conditions ---\n', opt.nTrials);
fprintf('  no control %6.2f (success %.2f)\n', b0.meanCens, b0.success);
fprintf('  ECR-II     %6.2f (success %.2f, hold %.2f)\n', b2.meanCens, b2.success, b2.hold);
fprintf('  ECR-III    %6.2f (success %.2f, hold %.2f)\n', b3.meanCens, b3.success, b3.hold);

resdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
if exist(resdir,'dir') ~= 7, mkdir(resdir); end
ecr_plot_clusters(m3, D, [1 2], fullfile(resdir,'lorenz_delay_clusters.png'));
ecr_save_model(fullfile(resdir,'lorenz_delay_models.mat'), m2, m3);
fprintf('\n=== done ===\n');
