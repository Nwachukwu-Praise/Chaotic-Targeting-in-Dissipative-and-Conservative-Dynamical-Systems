function audit = so_source_proxy_audit(cfg, route, sourceState)
%SO_SOURCE_PROXY_AUDIT Distances from the source/control segment to proxies.
if nargin < 3 || isempty(sourceState)
    sourceState = so_rectangle_center(cfg.sourceRectangle);
end
sourceState = sourceState(:);
sourceState(1) = so_wrap_x(sourceState(1));
samples = linspace(-cfg.controlAmplitude, cfg.controlAmplitude, max(101, cfg.curve.initialControlSamples));

rows = {};
minSourceDistance = Inf;
minControlledDistance = Inf;
zeroTimePossible = false;
for i = 1:numel(route.targetComponents)
    comps = route.targetComponents{i};
    for c = 1:numel(comps)
        comp = comps{c};
        sourceDistance = so_cylinder_distance(sourceState, comp.center);
        controlledDistances = zeros(numel(samples), 1);
        for s = 1:numel(samples)
            controlledDistances(s) = so_cylinder_distance(sourceState + [0; samples(s)], comp.center);
        end
        controlledDistance = min(controlledDistances);
        directAtZero = comp.contains(sourceState, cfg.containmentTolerance);
        directWithControl = any(comp.contains(sourceState + [zeros(1, numel(samples)); samples], ...
            cfg.containmentTolerance));
        minSourceDistance = min(minSourceDistance, sourceDistance);
        minControlledDistance = min(minControlledDistance, controlledDistance);
        zeroTimePossible = zeroTimePossible || directAtZero || directWithControl;
        rows(end + 1, :) = {i, string(comp.chainID), string(comp.id), comp.phasePointIndex, ...
            comp.center(1), comp.center(2), sourceDistance, controlledDistance, ...
            comp.radius, directAtZero, directWithControl}; %#ok<AGROW>
    end
end

names = {'routeIndex','chainID','componentID','phasePointIndex','proxyX','proxyY', ...
    'sourceDistance','controlledSegmentDistance','proxyRadius', ...
    'directAtZeroWithoutControl','directAtZeroWithAdmissibleControl'};
if isempty(rows)
    audit.byProxy = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    audit.byProxy = cell2table(rows, 'VariableNames', names);
end
audit.sourceState = sourceState;
audit.minimumSourceDistance = minSourceDistance;
audit.minimumControlledSegmentDistance = minControlledDistance;
audit.directContainmentPossibleAtN0 = zeroTimePossible;
audit.proxyRadius = cfg.proxyTargetRadius;
audit.sourceOutsideAllInitialProxies = minSourceDistance > cfg.proxyTargetRadius + cfg.containmentTolerance;
audit.controlSegmentOutsideAllInitialProxies = ...
    minControlledDistance > cfg.proxyTargetRadius + cfg.containmentTolerance;
end
