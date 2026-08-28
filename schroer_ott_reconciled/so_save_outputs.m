function paths = so_save_outputs(result)
%SO_SAVE_OUTPUTS Write the result MAT and the diagnostic CSV files.
cfg = result.configuration;
if ~exist(cfg.outputDirectory, 'dir')
    mkdir(cfg.outputDirectory);
end
paths = strings(0, 1);

name = char(cfg.caseName);
if isempty(name)
    name = 'case';
end

save(fullfile(cfg.outputDirectory, [name, '_result.mat']), 'result', '-v7.3');
paths(end + 1, 1) = string(fullfile(cfg.outputDirectory, [name, '_result.mat']));

paths = [paths; write(result.orbitCatalogue.byOmega, cfg, 'orbit_catalogue_by_omega.csv')];
paths = [paths; write(result.orbitCatalogue.byAbsResidue, cfg, 'orbit_catalogue_by_abs_residue.csv')];
paths = [paths; write(result.orbitCatalogue.diagnosticFractions, cfg, 'diagnostic_fractions.csv')];
paths = [paths; write(result.executedControls, cfg, 'executed_controls.csv')];
paths = [paths; write(result.equationOneDiagnostics, cfg, 'equation_one_diagnostics.csv')];
paths = [paths; write(so_stage_summary_table(result), cfg, 'stage_summary.csv')];
paths = [paths; write(execution_segment_table(result), cfg, 'execution_segments.csv')];
paths = [paths; write(result.resolutionFailures, cfg, 'resolution_failures.csv')];

if isfield(result, 'diagnosticManifolds') && ~isempty(result.diagnosticManifolds)
    paths = [paths; write(manifold_summary_table(result.diagnosticManifolds), cfg, ...
        'diagnostic_manifolds.csv')];
end
end

function p = write(tbl, cfg, fileName)
p = strings(0, 1);
if isempty(tbl) || width(tbl) == 0
    return;
end
target = fullfile(cfg.outputDirectory, fileName);
writetable(tbl, target);
p = string(target);
end

function tbl = execution_segment_table(result)
%EXECUTION_SEGMENT_TABLE One row per applied control.
%
% This is what the fundamental build could not report: where the state was
% immediately before and after each y-kick, which is the quantity a reader
% needs to check that |delta_n| never exceeded the bound and that the kick
% really is instantaneous and purely in y.
rows = {};
for i = 1:numel(result.executionSegments)
    s = result.executionSegments(i);
    pre = so_to_cylinder(s.preControlState);
    post = so_to_cylinder(s.postControlState);
    fin = so_to_cylinder(s.finalState);
    rows(end + 1, :) = {s.stage, s.targetID, s.control, s.iterations, ...
        pre(1), pre(2), post(1), post(2), fin(1), fin(2), ...
        s.consistencyError, s.isFinal}; %#ok<AGROW>
end
names = {'stage','targetID','controlY','iterations','preControlX','preControlY', ...
    'postControlX','postControlY','segmentEndX','segmentEndY', ...
    'propagationConsistencyError','isFinalSegment'};
if isempty(rows)
    tbl = cell2table(cell(0, numel(names)), 'VariableNames', names);
else
    tbl = cell2table(rows, 'VariableNames', names);
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
