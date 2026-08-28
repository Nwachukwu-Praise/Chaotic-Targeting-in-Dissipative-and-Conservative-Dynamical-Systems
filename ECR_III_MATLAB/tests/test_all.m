function nfail = test_all(quick)
%TEST_ALL  Run the whole ECR test suite.
%
%   TEST_ALL        full suite (includes a Lorenz pipeline, a few minutes)
%   TEST_ALL(true)  quick suite (skips the Lorenz pipeline test)
%
%   Returns the number of failed checks; prints a summary.  Run
%   STARTUP_ECR first.

if nargin < 1 || isempty(quick), quick = false; end
more off

tests = {@test_in_pi, @test_rbf, @test_cluster, @test_nmd, @test_control, ...
         @test_lorenz_system};
if ~quick
    tests{end+1} = @test_pipeline_logistic;
    tests{end+1} = @test_lorenz_pipeline;
else
    tests{end+1} = @test_pipeline_logistic;
end

nfail = 0;
for k = 1:numel(tests)
    name = func2str(tests{k});
    t0 = tic;
    try
        f = tests{k}();
        nfail = nfail + f;
        if f == 0
            fprintf('PASS  %-26s (%.1f s)\n', name, toc(t0));
        else
            fprintf('FAIL  %-26s (%d checks failed)\n', name, f);
        end
    catch err
        nfail = nfail + 1;
        fprintf('ERROR %-26s : %s\n', name, err.message);
    end
end
fprintf('\n%d failure(s)\n', nfail);
if nargout == 0, clear nfail, end
end
