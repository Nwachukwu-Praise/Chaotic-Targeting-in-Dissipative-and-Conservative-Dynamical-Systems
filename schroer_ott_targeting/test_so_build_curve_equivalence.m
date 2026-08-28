function results = test_so_build_curve_equivalence()
%TEST_SO_BUILD_CURVE_EQUIVALENCE Optimised builder versus verbatim reference.
%
%   The optimised so_build_curve changes only the evaluation strategy
%   (vectorised batch evaluation plus memoisation of curve points).  This
%   test asserts that, for open control segments and closed circular
%   boundaries, in both time directions and over a range of iterate counts,
%   the optimised builder reproduces the reference builder exactly:
%
%     * identical parameter sets;
%     * identical resolution status;
%     * identical subdivision depth and point count;
%     * identical curve points;
%     * identical reported maximum gap and midpoint deviation.
%
%   It also reports the speed-up, which is the reason the change was made.

cfg = schroer_ott_default_config();
cfg.verbose = false;

components = build_test_components(cfg);

rows = {};
for c = 1:numel(components)
    comp = components{c};
    for direction = [1, -1]
        for iterateCount = [0, 1, 3, 5, 6]
            tRef = tic;
            ref = so_build_curve_reference(comp.component, direction, iterateCount, cfg);
            refTime = toc(tRef);

            tNew = tic;
            [new, stats] = so_build_curve(comp.component, direction, iterateCount, cfg);
            newTime = toc(tNew);

            checks = compare_curves(ref, new);
            rows(end + 1, :) = { ...
                string(comp.name), direction, iterateCount, ...
                checks.parametersIdentical, checks.pointsIdentical, ...
                checks.statusIdentical, checks.depthIdentical, ...
                checks.gapIdentical, checks.midIdentical, ...
                checks.maxParameterDifference, checks.maxPointDifference, ...
                ref.pointCount, refTime, newTime, refTime / max(newTime, eps), ...
                stats.curvePointEvaluations, stats.cacheHits}; %#ok<AGROW>
        end
    end
end

results = cell2table(rows, 'VariableNames', {'component','direction','iterateCount', ...
    'parametersIdentical','pointsIdentical','statusIdentical','depthIdentical', ...
    'gapIdentical','midIdentical','maxParameterDifference','maxPointDifference', ...
    'pointCount','referenceSeconds','optimisedSeconds','speedup', ...
    'curvePointEvaluations','cacheHits'});

allIdentical = all(results.parametersIdentical) && all(results.pointsIdentical) && ...
    all(results.statusIdentical) && all(results.depthIdentical) && ...
    all(results.gapIdentical) && all(results.midIdentical);

fprintf('\n=== so_build_curve equivalence ===\n');
disp(results);
fprintf('cases                 : %d\n', height(results));
fprintf('all identical         : %d\n', allIdentical);
fprintf('max parameter diff    : %.3e\n', max(results.maxParameterDifference));
fprintf('max point diff        : %.3e\n', max(results.maxPointDifference));
fprintf('total reference time  : %.3f s\n', sum(results.referenceSeconds));
fprintf('total optimised time  : %.3f s\n', sum(results.optimisedSeconds));
fprintf('aggregate speedup     : %.1fx\n', ...
    sum(results.referenceSeconds) / max(sum(results.optimisedSeconds), eps));

if ~allIdentical
    error('SchroerOtt:BuildCurveEquivalence', ...
        'The optimised so_build_curve does not reproduce the reference builder.');
end
end

% -------------------------------------------------------------------------

function components = build_test_components(cfg)
components = {};

zSource = [0.884409; 0.506716];
components{end + 1} = struct('name', 'control-segment', ...
    'component', so_make_control_component(zSource, cfg));

catalogue = so_enumerate_periodic_orbits(cfg);
chains = catalogue.chains;
mask = strcmp({chains.classification}, 'direct-hyperbolic');
hyperbolic = chains(mask);
if isempty(hyperbolic)
    error('SchroerOtt:BuildCurveEquivalence', 'No direct-hyperbolic chain available.');
end
chain = hyperbolic(1);
components{end + 1} = struct('name', 'proxy-circle', ...
    'component', so_make_circle_component(chain.pointsLifted(:, 1), cfg.proxyTargetRadius, ...
    'equivalence-circle', chain.id, 1));

components{end + 1} = struct('name', 'final-rectangle', ...
    'component', so_make_rectangle_component(cfg.targetRectangle, cfg, 'equivalence-rectangle'));
end

function checks = compare_curves(ref, new)
pRef = ref.parameters(:).';
pNew = new.parameters(:).';
checks.parametersIdentical = numel(pRef) == numel(pNew) && isequal(pRef, pNew);
if numel(pRef) == numel(pNew)
    checks.maxParameterDifference = max([0, abs(pRef - pNew)]);
else
    checks.maxParameterDifference = Inf;
end

if isequal(size(ref.pointsLifted), size(new.pointsLifted))
    checks.pointsIdentical = isequal(ref.pointsLifted, new.pointsLifted);
    checks.maxPointDifference = max([0, abs(ref.pointsLifted(:) - new.pointsLifted(:)).']);
else
    checks.pointsIdentical = false;
    checks.maxPointDifference = Inf;
end

checks.statusIdentical = strcmp(string(ref.resolutionStatus), string(new.resolutionStatus));
checks.depthIdentical = ref.subdivisionDepth == new.subdivisionDepth && ...
    ref.pointCount == new.pointCount;
checks.gapIdentical = ref.maximumGap == new.maximumGap;
checks.midIdentical = ref.maximumMidpointDeviation == new.maximumMidpointDeviation;
end
