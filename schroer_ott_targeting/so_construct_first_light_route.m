function route = so_construct_first_light_route(catalogue, cfg)
%SO_CONSTRUCT_FIRST_LIGHT_ROUTE Rotation-bracket route in selected y-lift.
sourceCenter = so_rectangle_center(cfg.sourceRectangle);
targetRect = cfg.targetRectangle;
targetRect.yMin = targetRect.yMin + cfg.transport.targetLiftShift;
targetRect.yMax = targetRect.yMax + cfg.transport.targetLiftShift;
targetCenter = so_rectangle_center(targetRect);
dy = sign(targetCenter(2) - sourceCenter(2));
if dy == 0
    dy = 1;
end

chains = catalogue.chains;
directMask = strcmp({chains.classification}, 'direct-hyperbolic');
chains = chains(directMask);
if isempty(chains)
    error('SchroerOtt:NoDirectHyperbolicChains', ...
        'No direct-hyperbolic chains were enumerated.');
end
omega = [chains.omega];
lo = min(sourceCenter(2), targetCenter(2));
hi = max(sourceCenter(2), targetCenter(2));
margin = get_route_margin(cfg);
mask = route_bracket_mask(omega, lo, hi, margin, cfg);
chains = chains(mask);
if isempty(chains)
    error('SchroerOtt:NoRouteChains', ...
        'No direct-hyperbolic chains lie inside the declared y-lift bracket.');
end
if dy > 0
    [~, idx] = sort([chains.omega], 'ascend');
else
    [~, idx] = sort([chains.omega], 'descend');
end
chains = chains(idx);
if numel(chains) > cfg.maxStages
    error('SchroerOtt:RouteExceedsStageCap', ...
        'The first-light route has %d chains, exceeding cfg.maxStages=%d.', ...
        numel(chains), cfg.maxStages);
end

route.description = 'first-light rotation-bracket route';
route.sourceCenter = sourceCenter;
route.targetCenter = targetCenter;
route.transportDirection = dy;
route.bracketLower = lo;
route.bracketUpper = hi;
route.bracketMargin = margin;
route.bracketLowerBoundary = string(get_route_boundary(cfg, 'lowerBoundary'));
route.bracketUpperBoundary = string(get_route_boundary(cfg, 'upperBoundary'));
route.chains = chains;
route.rotationNumbers = [chains.omega];
route.targetComponents = cell(1, numel(chains));
for i = 1:numel(chains)
    route.targetComponents{i} = so_components_for_chain(chains(i), cfg);
end
route.finalTargetComponents = {so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target')};
route.ambiguities = detect_ambiguities(chains);
end

function margin = get_route_margin(cfg)
if isfield(cfg, 'routeBracket') && isfield(cfg.routeBracket, 'margin')
    margin = cfg.routeBracket.margin;
else
    margin = 1e-10;
end
end

function boundary = get_route_boundary(cfg, fieldName)
if isfield(cfg, 'routeBracket') && isfield(cfg.routeBracket, fieldName)
    boundary = cfg.routeBracket.(fieldName);
else
    boundary = 'open';
end
end

function mask = route_bracket_mask(omega, lo, hi, margin, cfg)
lowerBoundary = get_route_boundary(cfg, 'lowerBoundary');
upperBoundary = get_route_boundary(cfg, 'upperBoundary');
switch lower(lowerBoundary)
    case 'closed'
        lowerMask = omega >= lo - margin;
    case 'open'
        lowerMask = omega > lo + margin;
    otherwise
        error('SchroerOtt:InvalidRouteBracketPolicy', ...
            'Unknown lower route-bracket policy: %s.', lowerBoundary);
end
switch lower(upperBoundary)
    case 'closed'
        upperMask = omega <= hi + margin;
    case 'open'
        upperMask = omega < hi - margin;
    otherwise
        error('SchroerOtt:InvalidRouteBracketPolicy', ...
            'Unknown upper route-bracket policy: %s.', upperBoundary);
end
mask = lowerMask & upperMask;
end

function comps = so_components_for_chain(chain, cfg)
comps = cell(1, chain.period);
for i = 1:chain.period
    comps{i} = so_make_circle_component(chain.pointsLifted(:, i), cfg.proxyTargetRadius, ...
        sprintf('%s-phase-%02d', chain.id, i), chain.id, i);
end
end

function tbl = detect_ambiguities(chains)
if isempty(chains)
    tbl = table();
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
    tbl = table();
else
    tbl = cell2table(rows, 'VariableNames', {'omega','directChainCount','chainIDs'});
end
end
