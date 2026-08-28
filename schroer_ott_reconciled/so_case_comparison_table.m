function tbl = so_case_comparison_table(results)
%SO_CASE_COMPARISON_TABLE One row per case study, side by side.
rows = {};
for i = 1:numel(results)
    r = results{i};
    cfg = r.configuration;
    src = so_rectangle_center(cfg.sourceRectangle);
    tgt = so_rectangle_center(cfg.targetRectangle);
    if isempty(r.route.rotationNumbers)
        omegaText = "(empty bracket)";
    else
        omegaText = string(mat2str(r.route.rotationNumbers, 5));
    end
    if r.numberOfControls > 0
        maxControl = max(abs(r.executedControls.controlY));
        totalControl = sum(abs(r.executedControls.controlY));
        maxConsistency = max(r.executedControls.propagationConsistencyError);
    else
        maxControl = 0; totalControl = 0; maxConsistency = 0;
    end
    resid = -Inf;
    certified = true;
    for s = 1:numel(r.stagePlans)
        c = r.stagePlans(s).provisionalConnection;
        if c.success
            resid = max(resid, c.intersectionResidual);
            certified = certified && c.timeMinimumCertified;
        else
            certified = false;
        end
    end
    if ~isfinite(resid), resid = NaN; end
    if isfield(r, 'failureDiagnostics') && ~isempty(r.failureDiagnostics)
        failure = string(r.failureDiagnostics(end).failureCategory);
    else
        failure = "";
    end
    rows(end + 1, :) = { ...
        string(cfg.caseName), string(cfg.caseLabel), ...
        src(1), src(2), tgt(1), tgt(2), ...
        omegaText, numel(r.route.rotationNumbers), ...
        r.numberOfControls, r.totalExecutedIterations, ...
        maxControl, totalControl, cfg.controlAmplitude, ...
        r.finalState(1), r.finalState(2), r.targetContained, r.targetReached, ...
        resid, maxConsistency, certified, height(r.resolutionFailures), ...
        failure, r.runtimeSeconds}; %#ok<AGROW>
end
names = {'caseName','caseLabel','sourceX','sourceY','targetX','targetY', ...
    'routeOmegas','proxyCount','numberOfControls','totalIterations', ...
    'maxAbsControl','totalAbsControl','controlBound','finalX','finalY', ...
    'finalContainment','targetReached','maxIntersectionResidual', ...
    'maxConsistencyError','allMinimaCertified','unresolvedSplits', ...
    'failureCategory','runtimeSeconds'};
if isempty(rows)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    tbl = cell2table(rows, 'VariableNames', names);
end
end
