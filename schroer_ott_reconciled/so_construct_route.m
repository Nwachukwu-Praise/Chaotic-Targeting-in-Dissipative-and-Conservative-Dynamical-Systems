function route = so_construct_route(catalogue, cfg)
%SO_CONSTRUCT_ROUTE Rotation-bracket route in the selected y-lift.
%
% This is the fundamental-folder route constructor with one behavioural
% change: an EMPTY rotation bracket is a legitimate route rather than an
% error.  It is needed for the horizontal case study, where the source and
% the target lie in the same transport band, so no resonance separates them
% and the pass-targeting method degenerates -- correctly -- to a single
% Shinbrot forward-backward step onto the final target.
%
% Set cfg.route.allowEmptyBracket = false to restore the original error.

sourceCenter = so_rectangle_center(cfg.sourceRectangle);
targetRect = cfg.targetRectangle;
targetRect.yMin = targetRect.yMin + cfg.transport.targetLiftShift;
targetRect.yMax = targetRect.yMax + cfg.transport.targetLiftShift;
targetCenter = so_rectangle_center(targetRect);
dy = sign(targetCenter(2) - sourceCenter(2));
if dy == 0
    dy = 1;
end

allowEmpty = true;
if isfield(cfg, 'route') && isfield(cfg.route, 'allowEmptyBracket')
    allowEmpty = cfg.route.allowEmptyBracket;
end

chains = catalogue.chains;
directMask = strcmp({chains.classification}, 'direct-hyperbolic');
chains = chains(directMask);
if isempty(chains)
    error('SchroerOtt:NoDirectHyperbolicChains', ...
        'No direct-hyperbolic chains were enumerated.');
end
allDirect = chains;
omega = [chains.omega];
lo = min(sourceCenter(2), targetCenter(2));
hi = max(sourceCenter(2), targetCenter(2));
mask = omega > lo & omega < hi;
chains = chains(mask);

bracketEmpty = isempty(chains);
if bracketEmpty && ~allowEmpty
    error('SchroerOtt:NoRouteChains', ...
        'No direct-hyperbolic chains lie inside the declared y-lift bracket.');
end

if ~bracketEmpty
    if dy > 0
        [~, idx] = sort([chains.omega], 'ascend');
    else
        [~, idx] = sort([chains.omega], 'descend');
    end
    chains = chains(idx);
    if numel(chains) > cfg.maxStages
        error('SchroerOtt:RouteExceedsStageCap', ...
            'The route has %d chains, exceeding cfg.maxStages=%d.', ...
            numel(chains), cfg.maxStages);
    end
end

route.description = 'rotation-bracket route';
route.sourceCenter = sourceCenter;
route.targetCenter = targetCenter;
route.transportDirection = dy;
route.bracket = [lo, hi];
route.bracketEmpty = bracketEmpty;
if bracketEmpty
    route.chains = repmat(allDirect(1), 0, 1);
    route.rotationNumbers = zeros(1, 0);
    route.targetComponents = cell(1, 0);
else
    route.chains = chains;
    route.rotationNumbers = [chains.omega];
    route.targetComponents = cell(1, numel(chains));
    for i = 1:numel(chains)
        route.targetComponents{i} = components_for_chain(chains(i), cfg);
    end
end
route.finalTargetComponents = {so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target')};
route.ambiguities = detect_ambiguities(route.chains);
end

function comps = components_for_chain(chain, cfg)
comps = cell(1, chain.period);
for i = 1:chain.period
    comps{i} = so_make_circle_component(chain.pointsLifted(:, i), cfg.proxyTargetRadius, ...
        sprintf('%s-phase-%02d', chain.id, i), chain.id, i);
end
end

function tbl = detect_ambiguities(chains)
if isempty(chains)
    tbl = cell2table(cell(0, 3), 'VariableNames', {'omega','directChainCount','chainIDs'});
    return;
end
omegaVals = unique([chains.omega]);
rows = {};
for i = 1:numel(omegaVals)
    mask = abs([chains.omega] - omegaVals(i)) < 1e-12;
    if sum(mask) > 1
        rows(end + 1, :) = {omegaVals(i), sum(mask), strjoin(string({chains(mask).id}), ',')}; %#ok<AGROW>
    end
end
if isempty(rows)
    tbl = cell2table(cell(0, 3), 'VariableNames', {'omega','directChainCount','chainIDs'});
else
    tbl = cell2table(rows, 'VariableNames', {'omega','directChainCount','chainIDs'});
end
end
