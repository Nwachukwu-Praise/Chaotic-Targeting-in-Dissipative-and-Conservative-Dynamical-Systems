function ecr_plot_clusters(model, D, dims, fname)
%ECR_PLOT_CLUSTERS  Show the control regions and their ECR-III clusters.
%
%   ECR_PLOT_CLUSTERS(MODEL,D)              plot in the first two state dims
%   ECR_PLOT_CLUSTERS(MODEL,D,[i j])        choose the plotted dimensions
%   ECR_PLOT_CLUSTERS(MODEL,D,DIMS,FNAME)   also save the figure to FNAME
%
%   Background: the attractor data of D.  Foreground: one colour per control
%   region T_i, with the 2-sigma ellipse of every cluster C_ij drawn from its
%   mean mu_ij and covariance - i.e. the "loci of the clusters ... in terms
%   of the cluster means and the normalised cluster variances" of Section 3.3.

if nargin < 3 || isempty(dims), dims = [1 2]; end
S = model.S;
if S.N < 2, dims = [1 1]; end

figure('visible','off');
hold on
plot(D.Z(dims(1),:), D.Z(dims(2),:), '.', 'Color', [0.82 0.82 0.85], ...
     'MarkerSize', 3);

cols = lines_fallback(numel(model.levels));
h = []; lbl = {};
for k = 1:numel(model.levels)
    L = model.levels(k);
    for j = 1:numel(L.CL)
        idx = L.CL(j).idx;
        mu  = L.CL(j).mu;
        E   = ellipse2(mu(dims), L.CL(j).Sigma(dims,dims), 2);
        p1  = plot(E(1,:), E(2,:), '-', 'Color', cols(k,:), 'LineWidth', 1.2);
        plot(mu(dims(1)), mu(dims(2)), 'x', 'Color', cols(k,:), 'MarkerSize', 8);
        if j == 1, h(end+1) = p1; lbl{end+1} = sprintf('T_%d', L.index); end %#ok<AGROW>
    end
end
plot(S.zstar(dims(1)), S.zstar(dims(2)), 'kp', 'MarkerSize', 12, ...
     'MarkerFaceColor', 'k');
th = linspace(0,2*pi,100);
plot(S.zstar(dims(1)) + S.delta*cos(th), S.zstar(dims(2)) + S.delta*sin(th), ...
     'k-', 'LineWidth', 1.5);

xlabel(sprintf('z_%d', dims(1)));
ylabel(sprintf('z_%d', dims(2)));
title(sprintf('%s control regions - %s', model.variant, S.name), 'Interpreter', 'none');
if ~isempty(h), legend(h, lbl, 'Location', 'best'); end
grid on
box on

if nargin >= 4 && ~isempty(fname)
    print(gcf, fname, '-dpng', '-r150');
end
end

% =========================================================================
function E = ellipse2(mu, Sig, k)
[V, Dg] = eig((Sig+Sig')/2);
d = sqrt(max(diag(Dg), 0));
t = linspace(0, 2*pi, 80);
E = repmat(mu(:),1,numel(t)) + k*V*[d(1)*cos(t); d(2)*sin(t)];
end

% =========================================================================
function C = lines_fallback(n)
base = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; ...
        0.49 0.18 0.56; 0.47 0.67 0.19; 0.30 0.75 0.93; 0.64 0.08 0.18];
C = base(mod(0:n-1, size(base,1))+1, :);
end
