function so_save_outputs(result)
%SO_SAVE_OUTPUTS Save result MAT and diagnostic CSV files.
cfg = result.configuration;
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
save(fullfile(cfg.outputDirectory, 'schroer_ott_first_light_result.mat'), 'result');
writetable(result.orbitCatalogue.byOmega, fullfile(cfg.outputDirectory, 'orbit_catalogue_by_omega.csv'));
writetable(result.orbitCatalogue.byAbsResidue, fullfile(cfg.outputDirectory, 'orbit_catalogue_by_abs_residue.csv'));
writetable(result.orbitCatalogue.diagnosticFractions, fullfile(cfg.outputDirectory, 'diagnostic_fractions.csv'));
writetable(result.executedControls, fullfile(cfg.outputDirectory, 'executed_controls.csv'));
writetable(result.equationOneDiagnostics, fullfile(cfg.outputDirectory, 'equation_one_diagnostics.csv'));
if ~isempty(result.resolutionFailures) && width(result.resolutionFailures) > 0
    writetable(result.resolutionFailures, fullfile(cfg.outputDirectory, 'resolution_failures.csv'));
end
stageSummary = so_stage_summary_table(result);
writetable(stageSummary, fullfile(cfg.outputDirectory, 'stage_summary.csv'));
if isfield(result, 'diagnosticManifolds')
    writetable(manifold_summary_table(result.diagnosticManifolds), ...
        fullfile(cfg.outputDirectory, 'diagnostic_manifolds.csv'));
end
end

function tbl = manifold_summary_table(manifolds)
rows = {};
for i = 1:numel(manifolds)
    rows(end + 1, :) = {manifolds(i).chainID, manifolds(i).omega, ...
        manifolds(i).phasePointIndex, manifolds(i).branchType, manifolds(i).sign, ...
        manifolds(i).resolutionStatus, manifolds(i).pointCount, ...
        manifolds(i).maximumGap, manifolds(i).maximumMidpointDeviation}; %#ok<AGROW>
end
names = {'chainID','omega','phasePointIndex','branchType','sign', ...
    'resolutionStatus','pointCount','maximumGap','maximumMidpointDeviation'};
if isempty(rows)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    tbl = cell2table(rows, 'VariableNames', names);
end
end
