function plot_noise_results(results, target, params, noise)
%PLOT_NOISE_RESULTS Visualize bisection-only sigma_noise sensitivity.

if isfield(results, 'completed') && ~results.completed
    error('plot_noise_results requires a complete full noise result, not a pilot checkpoint.');
end

summary = results.summaryTable;
trials = results.trials;
modeList = unique(summary.Mode, 'stable');
colors = lines(max(numel(modeList), 1));

%% 1. Success fraction versus sigma_noise
figure('Name', 'Bisection noise success fraction', 'Color', 'w');
hold on;
for m = 1:numel(modeList)
    rows = strcmp(summary.Mode, modeList{m});
    y = summary.SuccessFraction(rows);
    lower = y - summary.WilsonLower95(rows);
    upper = summary.WilsonUpper95(rows) - y;
    errorbar(summary.SigmaNoise(rows), y, lower, upper, 'o-', ...
        'LineWidth', 1.5, 'MarkerFaceColor', colors(m, :), ...
        'Color', colors(m, :), 'DisplayName', modeList{m});
end
grid on;
xlabel('\sigma_{noise}: per-step, per-coordinate standard deviation');
ylabel('success fraction');
ylim([-0.05, 1.05]);
title('Bisection targeting reliability under Gaussian coordinate noise');
legend('Location', 'best');
hold off;

%% 2. Crossings to target for successful trials
figure('Name', 'Bisection noise crossings', 'Color', 'w');
hold on;
for m = 1:numel(modeList)
    mode = modeList{m};
    modeRows = strcmp(summary.Mode, mode);
    modeTrials = find(strcmp({trials.mode}, mode));

    for idx = modeTrials
        if trials(idx).hit
            plot(trials(idx).sigmaNoise, ...
                trials(idx).numCrossingsToTarget, 'o', ...
                'Color', colors(m, :), 'MarkerFaceColor', colors(m, :), ...
                'HandleVisibility', 'off');
        end
    end

    plot(summary.SigmaNoise(modeRows), ...
        summary.MeanCrossingsToTargetAmongSuccesses(modeRows), 's-', ...
        'LineWidth', 1.5, 'Color', colors(m, :), ...
        'MarkerFaceColor', colors(m, :), 'DisplayName', mode);
end
grid on;
xlabel('\sigma_{noise}: per-step, per-coordinate standard deviation');
ylabel('accepted crossings to target');
title('Successful noisy bisection-retargeted trials');
legend('Location', 'best');
hold off;

%% 3. Retargeting activity versus sigma_noise
figure('Name', 'Bisection noise retargeting activity', 'Color', 'w');
hold on;
for m = 1:numel(modeList)
    mode = modeList{m};
    modeRows = strcmp(summary.Mode, mode);
    modeTrials = find(strcmp({trials.mode}, mode));
    for idx = modeTrials
        plot(trials(idx).sigmaNoise, ...
            trials(idx).numRetargetingEvents, 'o', ...
            'Color', colors(m, :), 'MarkerFaceColor', colors(m, :), ...
            'HandleVisibility', 'off');
    end
    plot(summary.SigmaNoise(modeRows), ...
        summary.MeanRetargetingEvents(modeRows), 's-', ...
        'LineWidth', 1.5, 'Color', colors(m, :), ...
        'MarkerFaceColor', colors(m, :), 'DisplayName', mode);
end
grid on;
xlabel('\sigma_{noise}: per-step, per-coordinate standard deviation');
ylabel('retargeting events');
title('Cycle-cadence bisection retargeting activity');
legend('Location', 'best');
hold off;

%% 4. Final target error versus sigma_noise
figure('Name', 'Bisection noise final target error', 'Color', 'w');
hold on;
for m = 1:numel(modeList)
    mode = modeList{m};
    modeRows = strcmp(summary.Mode, mode);
    modeTrials = find(strcmp({trials.mode}, mode));
    for idx = modeTrials
        plot(trials(idx).sigmaNoise, trials(idx).finalDistance, 'o', ...
            'Color', colors(m, :), 'MarkerFaceColor', colors(m, :), ...
            'HandleVisibility', 'off');
    end
    plot(summary.SigmaNoise(modeRows), ...
        summary.MedianFinalTargetError(modeRows), 's-', ...
        'LineWidth', 1.5, 'Color', colors(m, :), ...
        'MarkerFaceColor', colors(m, :), 'DisplayName', mode);
end
yline(target.tolerance, 'k--', 'target tolerance');
grid on;
xlabel('\sigma_{noise}: per-step, per-coordinate standard deviation');
ylabel('|X - X_t| at final accepted crossing');
title('Noisy bisection final target error');
legend('Location', 'best');
hold off;

%% 5. Low-noise successful trajectory
lowIndex = find_low_noise_success(trials);
figure('Name', 'Low-noise bisection trajectory', 'Color', 'w');
plot_trial_trajectory(trials, lowIndex, target, params, ...
    'Low-noise successful bisection-controlled trajectory');

%% 6. High-noise failed or marginal trajectory
highIndex = find_high_noise_failed_or_marginal(trials);
figure('Name', 'High-noise bisection trajectory', 'Color', 'w');
plot_trial_trajectory(trials, highIndex, target, params, ...
    'High-noise failed or marginal bisection-controlled trajectory');

%% 7. Section sequences for representative noise levels
figure('Name', 'Bisection noisy section sequences', 'Color', 'w');
hold on;
selected = representative_trial_indices( ...
    trials, results.sigmaNoiseValues, modeList);
seqColors = lines(max(numel(selected), 1));
for k = 1:numel(selected)
    idx = selected(k);
    xValues = trials(idx).xSequence;
    crossingNumbers = 0:(numel(xValues) - 1);
    plot(crossingNumbers, xValues, 'o-', ...
        'LineWidth', 1.1, 'Color', seqColors(k, :), ...
        'DisplayName', sprintf('%s, sigma_{noise} = %.3g', ...
        trials(idx).mode, trials(idx).sigmaNoise));
end
yline(target.x, 'k-', 'X_t', 'LineWidth', 1.2);
yline(target.x + target.tolerance, 'k--');
yline(target.x - target.tolerance, 'k--');
grid on;
xlabel('accepted crossing number n');
ylabel('section coordinate X_n');
title('Noisy bisection-controlled Poincare sequences');
legend('Location', 'best');
hold off;

%% 8. Selected pControl values under noisy retargeting
figure('Name', 'Bisection noise pControl values', 'Color', 'w');
hold on;
for m = 1:numel(modeList)
    modeTrials = find(strcmp({trials.mode}, modeList{m}));
    for idx = modeTrials
        pValues = trials(idx).selectedPValues;
        if isempty(pValues)
            continue;
        end
        xValues = trials(idx).sigmaNoise * ones(size(pValues));
        plot(xValues, pValues, '.', 'Color', colors(m, :), ...
            'HandleVisibility', 'off');
    end
    plot(NaN, NaN, '.', 'Color', colors(m, :), ...
        'DisplayName', modeList{m});
end
yline(results.control.deltaP, 'r--', '+\Delta p');
yline(-results.control.deltaP, 'r--', '-\Delta p');
grid on;
xlabel('\sigma_{noise}: per-step, per-coordinate standard deviation');
ylabel('selected pControl');
title('Bisection pControl values used during noisy retargeting');
legend('Location', 'best');
hold off;

fprintf(['Noise plots use per-step, per-coordinate Gaussian standard ', ...
    'deviation and fixed RK4 step dt = %.4g.\n'], noise.dt);
fprintf('Directly comparable to Shinbrot Fig. 6 axis: false\n');
end

function plot_trial_trajectory(trials, idx, target, params, titleText)
if isempty(idx) || idx < 1 || idx > numel(trials) || isempty(trials(idx).X)
    axis off;
    text(0.05, 0.5, 'No suitable example trial was available.', ...
        'Units', 'normalized');
    title(titleText);
    return;
end

X = trials(idx).X;
plot3(X(:, 1), X(:, 2), X(:, 3), '-', ...
    'LineWidth', 0.8, 'Color', [0.05 0.35 0.7]);
hold on;
plot3(target.state(1), target.state(2), target.state(3), ...
    'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
draw_section_plane(X, params);
grid on;
axis tight;
view(35, 24);
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('%s (%s, sigma_{noise} = %.3g, hit = %d)', ...
    titleText, trials(idx).mode, trials(idx).sigmaNoise, trials(idx).hit));
legend({'noisy controlled trajectory', 'target fixed point', ...
    'section z = 26.921'}, 'Location', 'best');
hold off;
end

function draw_section_plane(X, params)
if isempty(X)
    return;
end
xLimits = [min(X(:, 1)), max(X(:, 1))];
yLimits = [min(X(:, 2)), max(X(:, 2))];
if diff(xLimits) == 0
    xLimits = xLimits + [-1, 1];
end
if diff(yLimits) == 0
    yLimits = yLimits + [-1, 1];
end
[xx, yy] = meshgrid(linspace(xLimits(1), xLimits(2), 2), ...
    linspace(yLimits(1), yLimits(2), 2));
zz = params.zSection * ones(size(xx));
surf(xx, yy, zz, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
    'FaceColor', [0.7 0.7 0.7]);
end

function idx = find_low_noise_success(trials)
hasTrajectory = arrayfun(@(trial) ~isempty(trial.X), trials);
hitIndices = find([trials.hit] & hasTrajectory);
if isempty(hitIndices)
    hitIndices = find([trials.hit]);
    if isempty(hitIndices)
        idx = [];
        return;
    end
end
[~, order] = sort([trials(hitIndices).sigmaNoise]);
idx = hitIndices(order(1));
end

function idx = find_high_noise_failed_or_marginal(trials)
hasTrajectory = arrayfun(@(trial) ~isempty(trial.X), trials);
failed = find(~[trials.hit] & hasTrajectory);
if ~isempty(failed)
    [~, order] = sort([trials(failed).sigmaNoise], 'descend');
    idx = failed(order(1));
    return;
end
available = find(hasTrajectory);
if isempty(available)
    [~, idx] = max([trials.sigmaNoise]);
else
    [~, order] = sort([trials(available).sigmaNoise], 'descend');
    idx = available(order(1));
end
end

function selected = representative_trial_indices( ...
    trials, sigmaNoiseValues, modeList)
if isempty(trials)
    selected = [];
    return;
end
levelIndices = unique(round(linspace( ...
    1, numel(sigmaNoiseValues), min(3, numel(sigmaNoiseValues)))));
selected = [];
for m = 1:numel(modeList)
    for k = 1:numel(levelIndices)
        sigmaNoise = sigmaNoiseValues(levelIndices(k));
        candidates = find(strcmp({trials.mode}, modeList{m}) & ...
            abs([trials.sigmaNoise] - sigmaNoise) < 10*eps);
        if isempty(candidates)
            continue;
        end
        hitLocal = find([trials(candidates).hit], 1, 'first');
        if ~isempty(hitLocal)
            selected(end + 1) = candidates(hitLocal); %#ok<AGROW>
        else
            selected(end + 1) = candidates(1); %#ok<AGROW>
        end
    end
end
end
