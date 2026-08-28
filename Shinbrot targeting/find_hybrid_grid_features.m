function features = find_hybrid_grid_features(pGrid, residualValues, tolerance)
%FIND_HYBRID_GRID_FEATURES Retain all direct hits and adjacent brackets.

if numel(pGrid) ~= numel(residualValues)
    error('pGrid and residualValues must have equal lengths.');
end

pGrid = pGrid(:).';
residualValues = residualValues(:).';
finiteValues = isfinite(residualValues);
features.directHitIndices = find( ...
    finiteValues & abs(residualValues) <= tolerance);

finitePairs = finiteValues(1:end-1) & finiteValues(2:end);
lowerIndices = find(finitePairs & ...
    residualValues(1:end-1) .* residualValues(2:end) < 0);
features.bracketIndexPairs = [lowerIndices(:), lowerIndices(:) + 1];
features.bracketIntervals = [pGrid(lowerIndices(:)).', ...
    pGrid(lowerIndices(:) + 1).'];
end
