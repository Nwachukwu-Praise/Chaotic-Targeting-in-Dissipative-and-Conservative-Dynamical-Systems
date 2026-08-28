function route = so_construct_omega_route(catalogue, cfg, omegaTargets, description)
%SO_CONSTRUCT_OMEGA_ROUTE Build a route from declared rotation numbers.
if nargin < 4
    description = 'declared omega route';
end
chains = catalogue.chains;
selected = repmat(chains(1), 0, 1);
for i = 1:numel(omegaTargets)
    mask = abs([chains.omega] - omegaTargets(i)) < 1e-10 & ...
        strcmp({chains.classification}, 'direct-hyperbolic');
    if ~any(mask)
        error('SchroerOtt:OmegaRouteMissingChain', ...
            'No direct-hyperbolic chain was found for omega %.12g.', omegaTargets(i));
    end
    candidates = chains(mask);
    [~, idx] = min(abs([candidates.residue]));
    selected(end + 1) = candidates(idx); %#ok<AGROW>
end

sourceCenter = so_rectangle_center(cfg.sourceRectangle);
targetRect = cfg.targetRectangle;
targetRect.yMin = targetRect.yMin + cfg.transport.targetLiftShift;
targetRect.yMax = targetRect.yMax + cfg.transport.targetLiftShift;
targetCenter = so_rectangle_center(targetRect);
dy = sign(targetCenter(2) - sourceCenter(2));
if dy == 0
    dy = 1;
end
if dy < 0
    selected = fliplr(selected);
end

route.description = char(description);
route.sourceCenter = sourceCenter;
route.targetCenter = targetCenter;
route.transportDirection = dy;
route.bracketLower = min(sourceCenter(2), targetCenter(2));
route.bracketUpper = max(sourceCenter(2), targetCenter(2));
route.bracketMargin = cfg.routeBracket.margin;
route.bracketLowerBoundary = string(cfg.routeBracket.lowerBoundary);
route.bracketUpperBoundary = string(cfg.routeBracket.upperBoundary);
route.chains = selected;
route.rotationNumbers = [selected.omega];
route.targetComponents = cell(1, numel(selected));
for i = 1:numel(selected)
    route.targetComponents{i} = components_for_chain(selected(i), cfg);
end
route.finalTargetComponents = {so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target')};
route.ambiguities = table();
end

function comps = components_for_chain(chain, cfg)
comps = cell(1, chain.period);
for i = 1:chain.period
    comps{i} = so_make_circle_component(chain.pointsLifted(:, i), cfg.proxyTargetRadius, ...
        sprintf('%s-phase-%02d', chain.id, i), chain.id, i);
end
end
