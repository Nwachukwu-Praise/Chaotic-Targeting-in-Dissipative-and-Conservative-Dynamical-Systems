function results = test_schroer_ott_reconciled(varargin)
%TEST_SCHROER_OTT_RECONCILED Self-checks for the reconciled build.
%
%   results = TEST_SCHROER_OTT_RECONCILED()            % fast tests only
%   results = TEST_SCHROER_OTT_RECONCILED('Full', true)% also solves a case
%
% The fast set needs no targeting run and takes a few seconds.  The full set
% additionally solves the diagonal case and checks the executed schedule
% against an independent replay, which is the check that matters most: it
% re-derives the trajectory from the recorded controls alone, using nothing
% the solver stored, and compares endpoints.

opts = so_parse_options(varargin, struct('Full', false, 'Verbose', true));

tests = {};
tests{end + 1} = @test_map_inverse_roundtrip;
tests{end + 1} = @test_area_preservation;
tests{end + 1} = @test_wrap_helpers;
tests{end + 1} = @test_periodic_orbits;
tests{end + 1} = @test_paper_fractions_present;
tests{end + 1} = @test_rectangle_component;
tests{end + 1} = @test_circle_component;
tests{end + 1} = @test_curve_matches_direct_iteration;
tests{end + 1} = @test_intersection_backends_agree;
tests{end + 1} = @test_cylinder_wrappers;
tests{end + 1} = @test_route_brackets;
tests{end + 1} = @test_case_centres_are_chaotic;
tests{end + 1} = @test_phase_portrait_window;
tests{end + 1} = @test_interquartile_range;
tests{end + 1} = @test_uncontrolled_transport_trivial;
if opts.Full
    tests{end + 1} = @test_diagonal_case_and_replay;
    tests{end + 1} = @test_zero_noise_replay_is_exact;
    tests{end + 1} = @test_horizontal_route_is_single_stage;
end

rows = {};
for i = 1:numel(tests)
    name = func2str(tests{i});
    timer = tic;
    try
        detail = tests{i}();
        passed = true;
        message = string(detail);
    catch err
        passed = false;
        message = string(err.message);
    end
    rows(end + 1, :) = {i, string(name), passed, message, toc(timer)}; %#ok<AGROW>
    if opts.Verbose
        if passed
            fprintf('  [pass] %-42s %s\n', name, message);
        else
            fprintf('  [FAIL] %-42s %s\n', name, message);
        end
    end
end
results = cell2table(rows, 'VariableNames', {'index','test','passed','message','seconds'});
if opts.Verbose
    fprintf('\n%d of %d tests passed.\n', sum(results.passed), height(results));
end
end

% ======================================================================= %

function msg = test_map_inverse_roundtrip()
cfg = so_reconciled_config();
rng(7);
z = [rand(1, 200); 2 * rand(1, 200) - 0.5];
err = max(vecnorm(so_standard_map_inverse_lifted(so_standard_map_lifted(z, cfg), cfg) - z));
assert(err < 1e-12, 'round-trip error %g', err);
msg = sprintf('max round-trip error %.2e', err);
end

function msg = test_area_preservation()
cfg = so_reconciled_config();
rng(11);
worst = 0;
for i = 1:100
    z = [rand; rand];
    worst = max(worst, abs(det(so_jacobian(z, cfg)) - 1));
end
assert(worst < 1e-13, 'det J deviates by %g', worst);
msg = sprintf('max |det J - 1| = %.2e', worst);
end

function msg = test_wrap_helpers()
assert(abs(so_wrap_x(1.25) - 0.25) < 1e-15, 'so_wrap_x');
assert(abs(so_wrap_x(-0.25) - 0.75) < 1e-15, 'so_wrap_x negative');
assert(abs(so_wrap_diff_x(0.9) + 0.1) < 1e-15, 'so_wrap_diff_x');
d = so_cylinder_distance([0.02; 1], [0.98; 1]);
assert(abs(d - 0.04) < 1e-15, 'cylinder distance across the seam: %g', d);
msg = 'wrap, lift and cylinder-distance helpers agree';
end

function msg = test_periodic_orbits()
cfg = so_reconciled_config();
cat = so_enumerate_periodic_orbits(cfg);
assert(~isempty(cat.chains), 'no chains enumerated');
worst = 0;
for i = 1:numel(cat.chains)
    c = cat.chains(i);
    r = so_iterate(c.root, c.period, cfg, 1, false) - c.root - [c.winding; 0];
    worst = max(worst, norm(r));
    % every phase point must return to itself after one full period
    z = c.pointsLifted(:, 1);
    zEnd = so_iterate(z, c.period, cfg, 1, false) - [c.winding; 0];
    worst = max(worst, norm(zEnd - z));
end
assert(worst < 1e-9, 'worst periodic residual %g', worst);
msg = sprintf('%d chains, worst residual %.2e', numel(cat.chains), worst);
end

function msg = test_paper_fractions_present()
cfg = so_reconciled_config();
cat = so_enumerate_periodic_orbits(cfg);
t = cat.diagnosticFractions;
assert(all(t.represented), 'a Figure-3 rotation number is missing');
assert(all(t.directHyperbolic), 'a Figure-3 resonance has no direct-hyperbolic chain');
msg = sprintf('all %d Figure-3 fractions present and direct-hyperbolic', height(t));
end

function msg = test_rectangle_component()
cfg = so_reconciled_config();
comp = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
s = linspace(0, 1, 401);
pts = comp.gamma0(s);
assert(all(comp.contains(pts, 1e-12)), 'boundary points are not contained');
closure = norm(comp.gamma0(0) - comp.gamma0(1));
assert(closure < 1e-12, 'boundary does not close: %g', closure);
centre = so_rectangle_center(cfg.targetRectangle);
assert(comp.contains(centre, 0), 'centre not contained');
outside = centre + [0; 10];
assert(~comp.contains(outside, 0), 'a far point is reported contained');
msg = sprintf('rectangle boundary closes to %.2e and contains its samples', closure);
end

function msg = test_circle_component()
cfg = so_reconciled_config();
centre = [0.3; 0.62];
comp = so_make_circle_component(centre, cfg.proxyTargetRadius, 'c', 'chain-x', 1);
s = linspace(0, 1, 257);
pts = comp.gamma0(s);
r = so_cylinder_distance(pts, centre);
assert(max(abs(r - cfg.proxyTargetRadius)) < 1e-12, 'radius drift %g', ...
    max(abs(r - cfg.proxyTargetRadius)));
assert(comp.contains(centre, 0), 'centre not contained');
msg = sprintf('circle radius exact to %.2e', max(abs(r - cfg.proxyTargetRadius)));
end

function msg = test_curve_matches_direct_iteration()
% The adaptive curve must be the image of gamma0 under the same map used
% everywhere else; this recomputes a sample of it from scratch.
cfg = so_reconciled_config();
comp = so_make_control_component([0.665; 0.5], cfg);
n = 5;
curve = so_build_curve(comp, 1, n, cfg);
idx = round(linspace(1, numel(curve.parameters), min(25, numel(curve.parameters))));
worst = 0;
for i = idx
    s = curve.parameters(i);
    direct = so_iterate(comp.gamma0(s), n, cfg, 1, false);
    stored = curve.pointsLifted(:, i);
    worst = max(worst, abs(so_wrap_diff_x(direct(1) - stored(1))));
    worst = max(worst, abs(direct(2) - stored(2)));
end
assert(worst < 1e-10, 'curve deviates from direct iteration by %g', worst);
assert(strcmp(curve.resolutionStatus, 'resolved'), ...
    'curve unresolved at n=%d: %s', n, curve.resolutionStatus);
msg = sprintf('curve at n=%d matches direct iteration to %.2e (%d points)', ...
    n, worst, curve.pointCount);
end

function msg = test_intersection_backends_agree()
% so_find_polyline_intersections_indexed is the production path and may
% take a polyxpoly shortcut when the Mapping Toolbox is present.
% so_find_polyline_intersections_exhaustive is the brute-force reference.
% They must find the same crossings; this is the check that the spatial
% index (or the toolbox shortcut) never silently loses one.
cfg = so_reconciled_config();
src = so_make_control_component([0.665; 0.5], cfg);
tgt = so_make_circle_component([0.2326; 0.7326], cfg.proxyTargetRadius, ...
    'probe', 'chain-probe', 3);
% The property that matters is that the indexed path never LOSES a
% crossing the reference finds; extra near-duplicate hits at a shared
% vertex are harmless, so they are reported rather than failed.
compared = 0;
matched = 0;
extra = 0;
for nForward = 0:3
    fCurve = so_build_curve(src, 1, nForward, cfg);
    for nBackward = 0:6
        bCurve = so_build_curve(tgt, -1, nBackward, cfg);
        indexed = so_find_polyline_intersections_indexed(fCurve, bCurve, cfg);
        exhaustive = so_find_polyline_intersections_exhaustive(fCurve, bCurve, cfg);
        compared = compared + 1;
        if isempty(exhaustive)
            extra = extra + numel(indexed);
            continue;
        end
        assert(~isempty(indexed), ...
            'nF=%d nB=%d: reference found %d crossings, indexed found none', ...
            nForward, nBackward, numel(exhaustive));
        pIndexed = cell2mat({indexed.point});
        for h = 1:numel(exhaustive)
            p = exhaustive(h).point;
            d = min(hypot(pIndexed(1, :) - p(1), pIndexed(2, :) - p(2)));
            assert(d < 1e-6, ...
                'nF=%d nB=%d: reference crossing at (%.6g, %.6g) is %g from the nearest indexed one', ...
                nForward, nBackward, p(1), p(2), d);
            matched = matched + 1;
        end
        extra = extra + max(0, numel(indexed) - numel(exhaustive));
    end
end
msg = sprintf('%d splits compared, %d reference crossings all matched (%d extra)', ...
    compared, matched, extra);
end

function msg = test_cylinder_wrappers()
cfg = so_reconciled_config();
rng(13);
z = [rand(1, 50); rand(1, 50)];
a = so_standard_map_cylinder(z, cfg);
b = so_standard_map_lifted(z, cfg);
assert(max(abs(a(2, :) - b(2, :))) < 1e-15, 'cylinder wrapper altered y');
assert(max(abs(so_wrap_diff_x(a(1, :) - b(1, :)))) < 1e-12, ...
    'cylinder wrapper altered x beyond a whole period');
assert(all(a(1, :) >= 0 & a(1, :) < 1), 'cylinder wrapper left x outside [0,1)');
c = so_standard_map_inverse_cylinder(b, cfg);
assert(max(abs(so_wrap_diff_x(c(1, :) - z(1, :)))) < 1e-12, ...
    'inverse cylinder wrapper does not undo the map in x');
assert(max(abs(c(2, :) - z(2, :))) < 1e-12, ...
    'inverse cylinder wrapper does not undo the map in y');
msg = 'cylinder wrappers agree with the lifted map modulo one period in x';
end

function msg = test_route_brackets()
cfg = so_reconciled_config();
cat = so_enumerate_periodic_orbits(cfg);
expected = struct( ...
    'diagonal',   [3/5, 2/3], ...
    'horizontal', [], ...
    'vertical',   [1/2, 3/5, 2/3]);
names = fieldnames(expected);
for i = 1:numel(names)
    caseCfg = so_case_config(names{i});
    route = so_construct_route(cat, caseCfg);
    want = expected.(names{i});
    got = route.rotationNumbers;
    assert(numel(got) == numel(want), ...
        '%s: expected %d route resonances, got %d (%s)', names{i}, ...
        numel(want), numel(got), mat2str(got, 4));
    if ~isempty(want)
        assert(max(abs(sort(got) - sort(want))) < 1e-12, ...
            '%s: rotation numbers %s do not match %s', names{i}, ...
            mat2str(got, 6), mat2str(want, 6));
    end
    assert(route.bracketEmpty == isempty(want), ...
        '%s: bracketEmpty flag disagrees with the route', names{i});
end
msg = 'diagonal 2, horizontal 0, vertical 3 route resonances as designed';
end

function msg = test_case_centres_are_chaotic()
% Every source and target centre must lie in the connected chaotic region.
% A regular island would give a finite-time Lyapunov exponent near zero and
% would make the whole exercise meaningless, so this is worth asserting.
cfg = so_reconciled_config();
names = {'diagonal', 'horizontal', 'vertical'};
worstSmallest = Inf;
for i = 1:numel(names)
    c = so_case_config(names{i});
    for z0 = [so_rectangle_center(c.sourceRectangle), so_rectangle_center(c.targetRectangle)]
        lam = ftle(z0, cfg, 3000);
        assert(lam > 0.10, '%s: centre (%.4g, %.4g) has FTLE %.4g -- regular', ...
            names{i}, z0(1), z0(2), lam);
        worstSmallest = min(worstSmallest, lam);
    end
end
msg = sprintf('all six case centres chaotic, smallest FTLE %.3f per iterate', worstSmallest);
end

function lam = ftle(z0, cfg, n)
z = z0;
v = [1; 0];
s = 0;
for i = 1:n
    v = so_jacobian(z, cfg) * v;
    nv = norm(v);
    s = s + log(nv);
    v = v / nv;
    z = so_standard_map_lifted(z, cfg);
end
lam = s / n;
end

function msg = test_phase_portrait_window()
cfg = so_reconciled_config();
cfg.background.xSeedCount = 12;
cfg.background.ySeedCount = 8;
cfg.background.iterations = 40;
w = [0.3, 0.8];
pts = so_phase_portrait(cfg, w);
assert(~isempty(pts), 'phase portrait is empty');
assert(all(pts(1, :) >= 0 & pts(1, :) < 1), 'x not wrapped into [0,1)');
assert(all(pts(2, :) >= w(1) - 1e-12 & pts(2, :) <= w(2) + 1e-12), ...
    'points outside the requested y window');
msg = sprintf('%d background points, all inside the window', size(pts, 2));
end

function msg = test_interquartile_range()
v = so_interquartile_range([1 2 3 4]);
assert(abs(v - 2) < 1e-12, 'IQR of 1:4 should be 2, got %g', v);
assert(isnan(so_interquartile_range([])), 'empty input should give NaN');
assert(so_interquartile_range(5) == 0, 'single value should give 0');
msg = 'IQR helper matches the prctile convention';
end

function msg = test_uncontrolled_transport_trivial()
cfg = so_reconciled_config();
inside = so_rectangle_center(cfg.targetRectangle);
t = so_uncontrolled_transport_time(inside, cfg, 10);
assert(t == 0, 'a point already in the target should give 0, got %g', t);
far = [0.5; -8];
t2 = so_uncontrolled_transport_time(far, cfg, 25);
assert(isnan(t2), 'an unreachable point within the cap should give NaN, got %g', t2);
msg = 'uncontrolled transport time handles arrival and cap correctly';
end

function msg = test_diagonal_case_and_replay()
result = so_run_case('diagonal', 'Plot', false, 'Save', false, ...
    'Verbose', false, 'Manifolds', false);
assert(result.targetReached, 'diagonal case did not reach the target');
cfg = result.configuration;

% independent replay from the recorded controls only
z = so_rectangle_center(cfg.sourceRectangle);
z(1) = so_wrap_x(z(1));
total = 0;
for i = 1:numel(result.executionSegments)
    seg = result.executionSegments(i);
    assert(abs(seg.control) <= cfg.controlAmplitude + 1e-12, ...
        'control %g exceeds the bound %g', seg.control, cfg.controlAmplitude);
    z = z + [0; seg.control];
    z = so_iterate(z, seg.iterations, cfg, 1, false);
    total = total + seg.iterations;
end
err = so_cylinder_distance(z, result.finalState);
assert(err < 1e-10, 'independent replay differs from the stored endpoint by %g', err);
assert(total == result.totalExecutedIterations, ...
    'segment iterations sum to %d, result reports %d', total, result.totalExecutedIterations);

comp = so_make_rectangle_component(cfg.targetRectangle, cfg, 'final-target');
assert(comp.contains(z, cfg.containmentTolerance), 'replayed endpoint is not in the target');
msg = sprintf('diagonal reached in %d iterations, replay agrees to %.2e', total, err);
end

function msg = test_zero_noise_replay_is_exact()
result = so_run_case('diagonal', 'Plot', false, 'Save', false, ...
    'Verbose', false, 'Manifolds', false);
rep = so_replay_with_noise(result, 0, zeros(1, sum([result.executionSegments.iterations])));
assert(rep.targetContained, 'zero-noise replay left the target');
assert(rep.endpointDisplacement < 1e-12, ...
    'zero-noise replay displaced the endpoint by %g', rep.endpointDisplacement);
msg = sprintf('zero-noise replay reproduces the endpoint to %.2e', rep.endpointDisplacement);
end

function msg = test_horizontal_route_is_single_stage()
result = so_run_case('horizontal', 'Plot', false, 'Save', false, ...
    'Verbose', false, 'Manifolds', false);
assert(result.route.bracketEmpty, 'horizontal route should have an empty bracket');
assert(isempty(result.switchEvaluations), ...
    'an empty route should evaluate no switches, found %d', numel(result.switchEvaluations));
assert(numel(result.stagePlans) == 1, ...
    'an empty route should produce exactly one stage, found %d', numel(result.stagePlans));
assert(result.targetReached, 'horizontal case did not reach the target');
msg = sprintf('horizontal solved as one forward-backward step, %d iterations', ...
    result.totalExecutedIterations);
end
