function catalogue = so_enumerate_periodic_orbits(cfg)
%SO_ENUMERATE_PERIODIC_ORBITS Enumerate fixed-(p,m) lifted periodic chains.
chains = struct([]);
chainCount = 0;
for p = cfg.orbit.minPeriod:cfg.orbit.maxPeriod
    for m = 1:(p - 1)
        if gcd(p, m) ~= 1
            continue;
        end
        seeds = periodic_seeds(p, m, cfg);
        localChains = struct([]);
        for s = 1:size(seeds, 2)
            [root, residual, converged] = so_newton_periodic_orbit(seeds(:, s), p, m, cfg);
            if ~converged || residual > cfg.orbit.residualTolerance
                continue;
            end
            chain = so_periodic_chain_from_root(root, p, m, residual, cfg);
            if chain.isLowerPeriod
                continue;
            end
            if ~is_duplicate_chain(chain, localChains, cfg)
                localChains = [localChains, chain]; %#ok<AGROW>
            end
        end
        for j = 1:numel(localChains)
            chainCount = chainCount + 1;
            localChains(j).id = sprintf('chain-%03d', chainCount);
            chains = [chains, localChains(j)]; %#ok<AGROW>
        end
    end
end

for i = 1:numel(chains)
    keyMask = arrayfun(@(c) c.period == chains(i).period && c.winding == chains(i).winding, chains);
    chains(i).samePMCount = sum(keyMask);
end

catalogue.chains = chains;
catalogue.byOmega = chain_table(chains, 'omega');
catalogue.byAbsResidue = chain_table(chains, 'absResidue');
catalogue.diagnosticFractions = diagnostic_fraction_table(chains);
end

function seeds = periodic_seeds(p, m, cfg)
omega = m / p;
xSeeds = linspace(0, 1, cfg.orbit.seedXCount + 1);
xSeeds(end) = [];
ySeeds = omega + cfg.orbit.seedYOffsets;
% Include symmetry-line style seeds used frequently for the standard map.
xSeeds = unique([xSeeds, 0, 0.25, 0.5, 0.75]);
ySeeds = unique([ySeeds, omega, 0, 0.5, 1]);
seeds = zeros(2, numel(xSeeds) * numel(ySeeds));
idx = 0;
for x = xSeeds
    for y = ySeeds
        idx = idx + 1;
        seeds(:, idx) = [x; y];
    end
end
end

function tf = is_duplicate_chain(chain, chains, cfg)
tf = false;
for i = 1:numel(chains)
    if chain_distance(chain, chains(i)) < cfg.orbit.deduplicationTolerance
        tf = true;
        return;
    end
end
end

function d = chain_distance(a, b)
if a.period ~= b.period || a.winding ~= b.winding
    d = Inf;
    return;
end
A = a.pointsCylinder;
B = b.pointsCylinder;
p = a.period;
best = Inf;
for shift = 0:(p - 1)
    Bshift = B(:, mod((0:p - 1) + shift, p) + 1);
    dx = so_wrap_diff_x(A(1, :) - Bshift(1, :));
    dy = A(2, :) - Bshift(2, :);
    % Permit integer y-cell copies only for deduplication; keep lift in a.
    dy = dy - round(mean(dy));
    best = min(best, max(hypot(dx, dy)));
end
d = best;
end

function tbl = chain_table(chains, mode)
if isempty(chains)
    tbl = table();
    return;
end
id = string({chains.id}).';
period = [chains.period].';
winding = [chains.winding].';
omega = [chains.omega].';
trace = [chains.trace].';
residue = [chains.residue].';
absResidue = abs(residue);
residual = [chains.residual].';
averageY = [chains.averageY].';
averageYError = [chains.averageYError].';
classification = string({chains.classification}).';
samePMCount = [chains.samePMCount].';
tbl = table(id, period, winding, omega, trace, residue, absResidue, residual, ...
    averageY, averageYError, classification, samePMCount);
switch mode
    case 'omega'
        tbl = sortrows(tbl, {'omega','period','winding','trace'});
    case 'absResidue'
        tbl = sortrows(tbl, {'absResidue','omega'});
end
end

function tbl = diagnostic_fraction_table(chains)
fractions = [1/4; 1/3; 2/5; 1/2; 3/5; 2/3; 3/4];
represented = false(size(fractions));
directHyperbolic = false(size(fractions));
for i = 1:numel(fractions)
    mask = abs([chains.omega] - fractions(i)) < 1e-12;
    represented(i) = any(mask);
    if any(mask)
        directHyperbolic(i) = any(strcmp({chains(mask).classification}, 'direct-hyperbolic'));
    end
end
tbl = table(fractions, represented, directHyperbolic, ...
    'VariableNames', {'omega','represented','directHyperbolic'});
end
