function so_print_result_summary(result)
%SO_PRINT_RESULT_SUMMARY Console summary for one targeting case.
cfg = result.configuration;
src = cfg.sourceRectangle;
tgt = cfg.targetRectangle;
fprintf('\n---- deterministic targeting summary: %s ----\n', char(cfg.caseName));
fprintf('target reached: %d\n', result.targetReached);
fprintf('final containment: %d\n', result.targetContained);
fprintf('source rectangle: x=[%.6g, %.6g], y=[%.6g, %.6g]\n', ...
    src.xMin, src.xMax, src.yMin, src.yMax);
fprintf('target rectangle: x=[%.6g, %.6g], y=[%.6g, %.6g]\n', ...
    tgt.xMin, tgt.xMax, tgt.yMin, tgt.yMax);
fprintf('selected target lift: %d\n', cfg.transport.targetLiftShift);
fprintf('transport direction: %+d\n', result.route.transportDirection);
fprintf('k: %.12g\n', cfg.k);
fprintf('control amplitude: %.12g\n', cfg.controlAmplitude);
fprintf('iteration budget: nForward<=%d, nBackward<=%d, tau<=%d\n', ...
    cfg.maxForwardIterations, cfg.maxBackwardIterations, cfg.maxTotalTransferTime);
fprintf('maxPeriod: %d\n', cfg.orbit.maxPeriod);
if isempty(result.route.rotationNumbers)
    fprintf('route rotation numbers: (empty bracket, single forward-backward step)\n');
else
    fprintf('route rotation numbers: %s\n', mat2str(result.route.rotationNumbers, 6));
end
if isempty(result.skippedProxies)
    fprintf('skipped proxies: (none)\n');
else
    fprintf('skipped proxies: %s\n', strjoin(cellstr(result.skippedProxies), ', '));
end
fprintf('number of applied controls: %d\n', result.numberOfControls);
fprintf('executed map iterations: %d\n', result.totalExecutedIterations);
fprintf('final state: [%.12g, %.12g]\n', result.finalState(1), result.finalState(2));

if result.numberOfControls > 0
    maxControl = max(abs(result.executedControls.controlY));
    maxConsistency = max(result.executedControls.propagationConsistencyError);
else
    maxControl = 0;
    maxConsistency = 0;
end
resids = [];
for i = 1:numel(result.stagePlans)
    if result.stagePlans(i).provisionalConnection.success
        resids(end + 1) = result.stagePlans(i).provisionalConnection.intersectionResidual; %#ok<AGROW>
    end
end
fprintf('maximum applied abs(control): %.12g  (bound %.12g)\n', maxControl, cfg.controlAmplitude);
if isempty(resids)
    maxResidual = NaN;
else
    maxResidual = max(resids, [], 'omitnan');
end
fprintf('maximum intersection residual: %.12g\n', maxResidual);
fprintf('maximum propagation-consistency error: %.12g\n', maxConsistency);
fprintf('number of unresolved splits: %d\n', height(result.resolutionFailures));
cert = arrayfun(@(s) s.provisionalConnection.success && ...
    s.provisionalConnection.timeMinimumCertified, result.stagePlans);
selectedResolved = arrayfun(@(s) s.provisionalConnection.success, result.stagePlans);
fprintf('all selected connections resolved: %d\n', all(selectedResolved));
fprintf('all transfer-time minima certified: %d\n', all(cert));
fprintf('runtime seconds: %.3f\n', result.runtimeSeconds);

disp('Executed controls:');
disp(result.executedControls);

if ~result.targetReached && ~isempty(result.failureDiagnostics)
    fprintf('Failure category: %s (stage %d)\n', ...
        char(result.failureDiagnostics(end).failureCategory), ...
        result.failureDiagnostics(end).stage);
end
end
