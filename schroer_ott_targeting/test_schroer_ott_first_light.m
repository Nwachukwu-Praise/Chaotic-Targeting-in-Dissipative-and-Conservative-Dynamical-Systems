function testResults = test_schroer_ott_first_light()
%TEST_SCHROER_OTT_FIRST_LIGHT Focused Schroer-Ott validation tests.
cfg = schroer_ott_default_config();
cfg.saveFigures = false;
testResults = table(string.empty(0, 1), false(0, 1), string.empty(0, 1), ...
    'VariableNames', {'name','passed','message'});

catalogue = [];
route = [];
result = [];
constructedConnection = [];
constructedProfile = [];
constructedReplay = [];

tests = {
    'lifted forward/inverse consistency', @test_lifted_inverse
    'cylinder forward/inverse consistency', @test_cylinder_inverse
    'analytic Jacobian versus finite differences', @test_jacobian_fd
    'area preservation', @test_area
    'wrapped-x/lifted-x trajectory consistency', @test_wrapped_lifted_consistency
    'y remains unwrapped', @test_y_unwrapped
    'adaptive refinement of folded curve', @test_fold_refinement
    'midpoint deviation detects hidden fold', @test_midpoint_detection
    'closed-boundary wraparound refinement', @test_closed_wrap_refinement
    'indexed intersection versus exhaustive', @test_indexed_vs_exhaustive
    'explicit indexed intersection backend', @test_explicit_backend
    'crossing across x=0', @test_crossing_across_x_zero
    'known ordinary polyline crossing', @test_ordinary_crossing
    'positive-time constructed u0 recovery', @test_constructed_u0
    'zero-time containment minimises absolute control', @test_zero_time_min_abs
    'circular-target containment', @test_circle_containment
    'rectangle containment', @test_rectangle_containment
    'periodic-boundary rectangle containment', @test_periodic_rectangle
    'diagonal direct time ordering', @test_diagonal_direct
    'unresolved curve budget is recorded', @test_unresolved_budget_recorded
    'unresolved numerical status does not throw', @test_unresolved_no_abort
    'incremental family versus independent recomputation', @test_family_consistency
    'shared backward cache versus independent recomputation', @test_backward_cache
    'union target minimises over phase components', @test_union_component
    'phase-point tie-breaking', @test_phase_tie
    'diagonal search returns finite constructed result', @test_diagonal_constructed
    'fixed-(p,m) periodic residuals', @test_periodic_residuals
    'retained-chain average-y consistency', @test_average_y_consistency
    'multiple chains retained for same (p,m)', @test_multiple_chains
    'direct-hyperbolic omega=1/4 recovered', @test_direct_quarter
    'cyclic/cylinder deduplication', @test_deduplication
    'inverse-hyperbolic chains excluded from route', @test_inverse_excluded
    'elliptic chains excluded from route', @test_elliptic_excluded
    'near-parabolic classification', @test_near_parabolic
    'declared route bracket policy', @test_route_bracket_policy
    'source is outside all initial proxies', @test_source_proxy_audit
    'executed control bounds', @test_control_bounds
    'control bookkeeping columns are distinct', @test_control_bookkeeping
    'independent replay validates accepted stages', @test_replay_valid
    'replay detects corrupted planned path', @test_replay_corrupt_path
    'replay detects corrupted control', @test_replay_corrupt_control
    'j=0 switch mechanism is objective-based', @test_j_zero_mechanism
    'at least one positive switching segment executes', @test_positive_switch
    'exact pruning rule', @test_pruning_rule
    'executed-time accounting', @test_time_accounting
    'final target independent containment', @test_final_containment
    };

for testIdx = 1:size(tests, 1)
    name = string(tests{testIdx, 1});
    fn = tests{testIdx, 2};
    try
        fn();
        testResults = [testResults; table(name, true, "", ...
            'VariableNames', {'name','passed','message'})]; %#ok<AGROW>
    catch err
        testResults = [testResults; table(name, false, string(err.message), ...
            'VariableNames', {'name','passed','message'})]; %#ok<AGROW>
    end
end

disp(testResults);
if any(~testResults.passed)
    error('SchroerOtt:TestsFailed', '%d Schroer-Ott first-light tests failed.', sum(~testResults.passed));
end

    function test_lifted_inverse()
        z = [1.73; -0.42];
        z2 = so_standard_map_inverse_lifted(so_standard_map_lifted(z, cfg), cfg);
        assert(norm(z2 - z) < 1e-12);
    end

    function test_cylinder_inverse()
        z = [0.83; 1.27];
        z2 = so_standard_map_inverse_cylinder(so_standard_map_cylinder(z, cfg), cfg);
        assert(abs(so_wrap_diff_x(z2(1) - z(1))) < 1e-12);
        assert(abs(z2(2) - z(2)) < 1e-12);
    end

    function test_jacobian_fd()
        z = [0.31; 0.44];
        J = so_jacobian(z, cfg);
        h = 1e-7;
        e1 = [h; 0];
        e2 = [0; h];
        Jfd = [(so_standard_map_lifted(z + e1, cfg) - so_standard_map_lifted(z - e1, cfg)) / (2 * h), ...
            (so_standard_map_lifted(z + e2, cfg) - so_standard_map_lifted(z - e2, cfg)) / (2 * h)];
        assert(norm(J - Jfd) < 1e-7);
    end

    function test_area()
        assert(abs(det(so_jacobian([0.2; 0.3], cfg)) - 1) < 1e-12);
    end

    function test_wrapped_lifted_consistency()
        z = [2.2; 0.37];
        lifted = so_iterate(z, 7, cfg, 1, true);
        cyl = so_to_cylinder(lifted);
        zc = so_to_cylinder(z);
        for n = 1:7
            zc = so_standard_map_cylinder(zc, cfg);
        end
        assert(abs(so_wrap_diff_x(cyl(1, end) - zc(1))) < 1e-12);
        assert(abs(cyl(2, end) - zc(2)) < 1e-12);
    end

    function test_y_unwrapped()
        z = [0.1; 2.4];
        zn = so_standard_map_cylinder(z, cfg);
        assert(zn(2) > 2.0);
    end

    function test_fold_refinement()
        c = synthetic_component(@(s) [s; 0.05 * sin(8 * pi * s)], false, [0 1], 0.01);
        ccfg = cfg;
        ccfg.curve.initialControlSamples = 5;
        curve = so_build_curve(c, 1, 0, ccfg);
        assert(curve.pointCount > 5);
        assert(strcmp(curve.resolutionStatus, "resolved"));
    end

    function test_midpoint_detection()
        c = synthetic_component(@(s) [0.5 + 0.001 * cos(2 * pi * s); 0.5 + 0.08 * sin(2 * pi * s)], ...
            true, [0 1], 0.02);
        ccfg = cfg;
        ccfg.curve.initialBoundarySamples = 8;
        curve = so_build_curve(c, 1, 0, ccfg);
        assert(curve.pointCount > 8);
    end

    function test_closed_wrap_refinement()
        rect = struct('xMin', 0.98, 'xMax', 0.02, 'yMin', 0.1, 'yMax', 0.2, 'id', 'wrap-rect');
        comp = so_make_rectangle_component(rect, cfg, 'wrap-rect');
        curve = so_build_curve(comp, 1, 0, cfg);
        seg = so_curve_segments(curve);
        assert(size(seg.a, 2) == curve.pointCount);
        assert(strcmp(curve.resolutionStatus, "resolved"));
    end

    function test_indexed_vs_exhaustive()
        a = make_poly_curve([0.2 0.6; 0 1], [0 1], false);
        b = make_poly_curve([0.2 0.6; 1 0], [0 1], false);
        hi = so_find_polyline_intersections_indexed(a, b, cfg);
        he = so_find_polyline_intersections_exhaustive(a, b, cfg);
        assert(numel(hi) >= 1 && numel(he) >= 1);
        pts = [hi.point];
        assert(any(vecnorm(pts - [0.4; 0.5]) < 1e-12));
        assert(isfield(hi, 'sALo') && isfield(hi, 'sBHi'));
    end

    function test_explicit_backend()
        a = make_poly_curve([0.2 0.6; 0 1], [0 1], false);
        b = make_poly_curve([0.2 0.6; 1 0], [0 1], false);
        badCfg = cfg;
        badCfg.intersectionBackend = 'polyxpoly';
        didThrow = false;
        try
            so_find_polyline_intersections_indexed(a, b, badCfg);
        catch err
            didThrow = strcmp(err.identifier, 'SchroerOtt:IntersectionBackendMismatch');
        end
        assert(didThrow);
    end

    function test_crossing_across_x_zero()
        a = make_poly_curve([0.99 1.01; 0 1], [0 1], false);
        b = make_poly_curve([0 0; 1 0], [0 1], false);
        hi = so_find_polyline_intersections_indexed(a, b, cfg);
        assert(~isempty(hi));
    end

    function test_ordinary_crossing()
        h = so_segment_intersection([0;0], [1;1], [0;1], [1;0], 1e-12);
        assert(h.success && norm(h.point - [0.5; 0.5]) < 1e-12);
    end

    function test_constructed_u0()
        load_constructed_connection();
        assert(constructedConnection.success);
        assert(~constructedConnection.directContainment);
        assert(constructedConnection.totalIterations > 0);
        assert(constructedProfile.indexedCrossingQueries > 0);
        assert(constructedProfile.refinementCalls > 0);
        assert(abs(constructedConnection.control - 0.0007) < 5e-6);
        assert(constructedConnection.intersectionResidual < 1e-8);
        assert(constructedReplay.replayPassed);
        assert(constructedReplay.targetContained);
    end

    function test_zero_time_min_abs()
        z = [0.4; 0.5];
        c1Rect = struct('xMin', 0.39, 'xMax', 0.41, 'yMin', 0.4978, ...
            'yMax', 0.4980, 'id', 'large-negative');
        c2Rect = struct('xMin', 0.39, 'xMax', 0.41, 'yMin', 0.50035, ...
            'yMax', 0.50045, 'id', 'small-positive');
        c1 = so_make_rectangle_component(c1Rect, cfg, 'large-negative');
        c2 = so_make_rectangle_component(c2Rect, cfg, 'small-positive');
        [conn, ~] = so_resolve_connection(z, {c1, c2}, small_cfg(), so_new_performance_profile());
        assert(conn.success && conn.zeroTimePretest);
        assert(conn.targetComponent == "small-positive");
        assert(abs(conn.control - 0.00035) < 1e-12);
    end

    function test_circle_containment()
        comp = so_make_circle_component([0.99; 0.5], 0.03, 'circle', '', 1);
        assert(comp.contains([0.01; 0.5], 1e-12));
    end

    function test_rectangle_containment()
        comp = so_make_rectangle_component(cfg.targetRectangle, cfg, 'rect');
        assert(comp.contains([0.23; 0.72], 1e-12));
        assert(~comp.contains([0.5; 0.72], 1e-12));
    end

    function test_periodic_rectangle()
        rect = struct('xMin', 0.95, 'xMax', 0.05, 'yMin', 0, 'yMax', 1, 'id', 'wrap');
        assert(so_point_in_rectangle([0.99; 0.5], rect, 0));
        assert(so_point_in_rectangle([0.01; 0.5], rect, 0));
        assert(~so_point_in_rectangle([0.50; 0.5], rect, 0));
    end

    function test_diagonal_direct()
        z = [0.3; 0.3];
        comp = so_make_circle_component(z, 0.01, 'direct', '', 1);
        [conn, ~] = so_resolve_connection(z, {comp}, small_cfg(), so_new_performance_profile());
        assert(conn.success && conn.totalIterations == 0);
        assert(conn.zeroTimePretest);
    end

    function test_unresolved_budget_recorded()
        [conn, profile] = unresolved_budget_case();
        assert(~conn.success);
        assert(conn.resolutionHoles > 0);
        assert(profile.unresolvedSplits > 0);
        assert(~conn.timeMinimumCertified);
    end

    function test_unresolved_no_abort()
        [conn, ~] = unresolved_budget_case();
        assert(isstruct(conn));
        assert(isfield(conn, 'resolutionFailures'));
    end

    function test_family_consistency()
        comp = so_make_control_component([0.2; 0.3], cfg);
        [fam, ~] = so_build_curve_family(comp, 1, 2, small_cfg(), so_new_performance_profile(), 'forward');
        independent = so_build_curve(comp, 1, 2, small_cfg());
        assert(norm(fam{3}.pointsLifted(:, end) - independent.pointsLifted(:, end)) < 1e-12);
    end

    function test_backward_cache()
        comp = so_make_circle_component([0.4; 0.4], 0.02, 'cache', '', 1);
        ccfg = small_cfg();
        [cache, profile] = so_build_backward_cache({comp}, ccfg, so_new_performance_profile());
        [connA, ~] = so_resolve_connection([0.4; 0.4], {comp}, ccfg, profile, cache);
        [connB, ~] = so_resolve_connection([0.4; 0.4], {comp}, ccfg, so_new_performance_profile());
        assert(connA.success == connB.success);
    end

    function test_union_component()
        z = [0.5; 0.5];
        c1 = so_make_circle_component([0.7; 0.5], 0.01, 'far', '', 1);
        c2 = so_make_circle_component(z, 0.01, 'near', '', 2);
        [conn, ~] = so_resolve_connection(z, {c1, c2}, small_cfg(), so_new_performance_profile());
        assert(conn.success && conn.targetPhasePointIndex == 2);
    end

    function test_phase_tie()
        z = [0.5; 0.5];
        c1 = so_make_circle_component(z, 0.01, 'a', '', 1);
        c2 = so_make_circle_component(z, 0.01, 'b', '', 2);
        [conn, ~] = so_resolve_connection(z, {c1, c2}, small_cfg(), so_new_performance_profile());
        assert(conn.success && conn.targetPhasePointIndex == 1);
    end

    function test_diagonal_constructed()
        load_constructed_connection();
        assert(isfinite(constructedConnection.tauResolved));
        assert(constructedConnection.tauResolved == constructedConnection.totalIterations);
    end

    function test_periodic_residuals()
        load_catalogue();
        assert(max(catalogue.byOmega.residual) < cfg.orbit.residualTolerance * 1e3);
    end

    function test_average_y_consistency()
        load_catalogue();
        assert(max(catalogue.byOmega.averageYError) < 1e-10);
    end

    function test_multiple_chains()
        load_catalogue();
        mask = catalogue.byOmega.period == 4 & catalogue.byOmega.winding == 1;
        assert(sum(mask) >= 2);
    end

    function test_direct_quarter()
        load_catalogue();
        mask = abs(catalogue.byOmega.omega - 0.25) < 1e-12 & ...
            catalogue.byOmega.classification == "direct-hyperbolic";
        assert(any(mask));
    end

    function test_deduplication()
        load_catalogue();
        assert(numel(unique(catalogue.byOmega.id)) == height(catalogue.byOmega));
    end

    function test_inverse_excluded()
        load_catalogue();
        load_route();
        inverseIDs = catalogue.byOmega.id(catalogue.byOmega.classification == "inverse-hyperbolic");
        routeIDs = string({route.chains.id}).';
        assert(~isempty(inverseIDs));
        assert(~any(ismember(inverseIDs, routeIDs)));
    end

    function test_elliptic_excluded()
        load_catalogue();
        load_route();
        ellipticIDs = catalogue.byOmega.id(catalogue.byOmega.classification == "elliptic");
        routeIDs = string({route.chains.id}).';
        assert(~isempty(ellipticIDs));
        assert(~any(ismember(ellipticIDs, routeIDs)));
    end

    function test_near_parabolic()
        assert(strcmp(so_classify_trace(2 + 0.5 * cfg.orbit.traceTolerance, cfg.orbit.traceTolerance), ...
            'near-parabolic'));
    end

    function test_route_bracket_policy()
        load_route();
        assert(route.bracketLowerBoundary == "open");
        assert(route.bracketUpperBoundary == "open");
        assert(route.bracketMargin > 0);
        assert(all(route.rotationNumbers > route.bracketLower + route.bracketMargin));
        assert(all(route.rotationNumbers < route.bracketUpper - route.bracketMargin));
    end

    function test_source_proxy_audit()
        load_route();
        audit = so_source_proxy_audit(cfg, route);
        assert(audit.sourceOutsideAllInitialProxies);
        assert(audit.controlSegmentOutsideAllInitialProxies);
        assert(~audit.directContainmentPossibleAtN0);
    end

    function test_control_bounds()
        load_result();
        assert(max(abs(result.executedControls.controlY)) <= cfg.controlAmplitude + 10 * eps);
    end

    function test_control_bookkeeping()
        load_result();
        tbl = result.executedControls;
        required = ["preControlX","preControlY","postControlX","postControlY", ...
            "switchStateX","switchStateY","controlY","stageIterations","cumulativeIterations"];
        assert(all(ismember(required, string(tbl.Properties.VariableNames))));
        dy = tbl.postControlY - tbl.preControlY;
        assert(max(abs(dy - tbl.controlY)) < 1e-12);
        assert(any(abs(tbl.switchStateY - tbl.preControlY) > 1e-6));
    end

    function test_replay_valid()
        load_result();
        assert(all(result.executedControls.replayPassed));
        assert(max(result.executedControls.maxPathwiseReplayError) < cfg.propagationConsistencyTolerance);
    end

    function test_replay_corrupt_path()
        load_result();
        seg = first_positive_segment();
        planned = seg.replay.storedPlannedPath;
        planned(:, end) = planned(:, end) + [0; 1e-4];
        bad = so_replay_stage(seg.preControlState, seg.control, planned, seg.iterations, cfg, {});
        assert(~bad.replayPassed);
        assert(bad.maxPathwiseReplayError > cfg.propagationConsistencyTolerance);
    end

    function test_replay_corrupt_control()
        load_result();
        seg = first_positive_segment();
        bad = so_replay_stage(seg.preControlState, seg.control + 5e-5, ...
            seg.replay.storedPlannedPath, seg.iterations, cfg, {});
        assert(~bad.replayPassed);
        assert(bad.maxPathwiseReplayError > cfg.propagationConsistencyTolerance);
    end

    function test_j_zero_mechanism()
        load_result();
        selectedZero = false;
        for i = 1:numel(result.switchEvaluations)
            selected = result.switchEvaluations(i).selected;
            if selected.finite && selected.j == 0
                probes = result.switchEvaluations(i).probes;
                finite = [probes.finite];
                objectives = [probes.objectiveJ];
                jVals = [probes.j];
                assert(all(objectives(finite & jVals > 0) >= selected.objectiveJ));
                stage = result.stagePlans(i);
                assert(~stage.provisionalConnection.directContainment);
                selectedZero = true;
            end
        end
        assert(selectedZero);
    end

    function test_positive_switch()
        load_result();
        assert(any(result.executedControls.stageIterations > 0));
    end

    function test_pruning_rule()
        assert(~so_should_prune_probe(4, Inf));
        assert(~so_should_prune_probe(3, 4));
        assert(so_should_prune_probe(4, 4));
    end

    function test_time_accounting()
        load_result();
        assert(result.totalExecutedIterations == sum(result.executedControls.stageIterations));
    end

    function test_final_containment()
        load_result();
        finalComp = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
        finalReplay = result.executionSegments(end).replay;
        assert(result.targetContained);
        assert(finalReplay.targetContained);
        assert(finalComp.contains(result.finalState, cfg.containmentTolerance));
    end

    function load_catalogue()
        if isempty(catalogue)
            catalogue = so_enumerate_periodic_orbits(cfg);
        end
    end

    function load_route()
        load_catalogue();
        if isempty(route)
            route = so_construct_first_light_route(catalogue, cfg);
        end
    end

    function load_result()
        load_route();
        if isempty(result)
            result = so_multistage_targeting(cfg, catalogue, route);
        end
    end

    function load_constructed_connection()
        if ~isempty(constructedConnection)
            return;
        end
        ccfg = cfg;
        ccfg.maxForwardIterations = 4;
        ccfg.maxBackwardIterations = 4;
        ccfg.maxTotalTransferTime = 8;
        ccfg.curve.maxPoints = 5000;
        ccfg.curve.initialBoundarySamples = 32;
        z = [0.41; 0.27];
        u0 = 0.0007;
        zFinal = so_iterate(z + [0; u0], 4, ccfg, 1, false);
        comp = so_make_circle_component(zFinal, 1e-5, 'constructed', '', 7);
        constructedProfile = so_new_performance_profile();
        [constructedConnection, constructedProfile] = so_resolve_connection(z, {comp}, ...
            ccfg, constructedProfile);
        constructedReplay = so_replay_stage(z, constructedConnection.control, ...
            constructedConnection.plannedPath, constructedConnection.totalIterations, ccfg, {comp});
    end

    function [conn, profile] = unresolved_budget_case()
        ccfg = cfg;
        ccfg.maxForwardIterations = 2;
        ccfg.maxBackwardIterations = 2;
        ccfg.maxTotalTransferTime = 4;
        ccfg.curve.maxPoints = 4;
        ccfg.curve.initialBoundarySamples = 8;
        ccfg.curve.initialControlSamples = 5;
        comp = so_make_circle_component([0.73; 0.62], 0.02, 'budget', '', 1);
        profile = so_new_performance_profile();
        [conn, profile] = so_resolve_connection([0.2; 0.31], {comp}, ccfg, profile);
    end

    function seg = first_positive_segment()
        segments = result.executionSegments;
        idx = find(arrayfun(@(s) s.iterations > 0, segments), 1, 'first');
        assert(~isempty(idx));
        seg = segments(idx);
    end

    function ccfg = small_cfg()
        ccfg = cfg;
        ccfg.maxForwardIterations = 2;
        ccfg.maxBackwardIterations = 2;
        ccfg.maxTotalTransferTime = 4;
        ccfg.curve.maxPoints = 2000;
        ccfg.curve.initialBoundarySamples = 16;
    end
end

function component = synthetic_component(gamma, isClosed, range, charSize)
component.type = 'synthetic';
component.id = 'synthetic';
component.isClosed = isClosed;
component.parameterRange = range;
component.characteristicSize = charSize;
component.center = [NaN; NaN];
component.phasePointIndex = NaN;
component.chainID = "";
component.gamma0 = @(s) gamma(reshape(s, 1, []));
component.contains = @(z, tol) false(1, size(z, 2)); %#ok<INUSD>
end

function curve = make_poly_curve(points, parameters, isClosed)
curve.gamma0 = [];
curve.component = [];
curve.parameters = parameters;
curve.isClosed = isClosed;
curve.direction = 1;
curve.iterateCount = 0;
curve.pointsLifted = points;
curve.pointsCylinder = so_to_cylinder(points);
curve.resolutionStatus = "resolved";
curve.maximumGap = 0;
curve.maximumMidpointDeviation = 0;
curve.subdivisionDepth = 0;
curve.pointCount = size(points, 2);
curve.hMax = 0.01;
curve.hMid = 0.01;
end
