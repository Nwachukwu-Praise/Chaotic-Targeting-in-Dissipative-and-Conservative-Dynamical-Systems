function ecr_plot_run(S, R, fname)
%ECR_PLOT_RUN  Plot one closed-loop run produced by ECR_SIMULATE.
%
%   ECR_PLOT_RUN(S,R)         three stacked panels
%   ECR_PLOT_RUN(S,R,FNAME)   also save the figure as PNG
%
%   Panel 1: distance ||z_n - z*|| with the delta threshold, on a log scale.
%   Panel 2: applied parameter perturbations relative to their allowed range.
%   Panel 3: the control region that was activated at each step.

n = size(R.Z,2);
d = sqrt(sum((R.Z - repmat(S.zstar,1,n)).^2, 1));

figure('visible','off');

subplot(3,1,1);
semilogy(0:n-1, max(d,1e-12), 'o-', 'LineWidth', 1.1, 'MarkerSize', 4);
hold on
semilogy([0 n-1], S.delta*[1 1], 'r--', 'LineWidth', 1.2);
ylabel('||z_n - z^*||');
title(sprintf('%s: reaching time = %s steps', S.name, num2str(R.reach)), 'Interpreter', 'none');
grid on

subplot(3,1,2);
dp = (R.P - repmat(S.pnom,1,n)) ./ repmat(S.dpmax,1,n);
plot(0:n-1, dp', 'o-', 'LineWidth', 1.0, 'MarkerSize', 3);
hold on
plot([0 n-1], [1 1], 'k:', [0 n-1], [-1 -1], 'k:');
ylabel('\deltap / \deltap_{max}');
ylim([-1.3 1.3]);
legend(S.pnames, 'Location', 'best');
grid on

subplot(3,1,3);
stairs(0:n-1, R.reg, 'LineWidth', 1.2);
ylabel('active region');
xlabel('step n');
ylim([-1.5 max(2, max(R.reg)+0.5)]);
grid on

if nargin >= 3 && ~isempty(fname)
    print(gcf, fname, '-dpng', '-r150');
end
end
