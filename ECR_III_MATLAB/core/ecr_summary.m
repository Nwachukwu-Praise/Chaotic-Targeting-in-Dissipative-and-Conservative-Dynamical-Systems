function txt = ecr_summary(model)
%ECR_SUMMARY  Human readable description of a trained ECR model.
%
%   TXT = ECR_SUMMARY(MODEL) returns (and, if called without an output
%   argument, prints) a table of the control regions, their clusters, their
%   sizes and their training errors.

L = model.levels;
lines = {};
lines{end+1} = sprintf('%s model for "%s"', model.variant, model.S.name);
lines{end+1} = sprintf('  target z*      : [%s]', strtrim(sprintf('%.4f ', model.zstar)));
lines{end+1} = sprintf('  delta          : %g', model.delta);
lines{end+1} = sprintf('  regions        : %d (T0 ... T%d)', numel(L), L(end).index);
lines{end+1} = sprintf('  networks       : %d', model.info.nNets);
lines{end+1} = sprintf('  training time  : %.2f s', model.info.trainTime);
lines{end+1} = '  ------------------------------------------------------------';
lines{end+1} = '   region   points   radius r   clusters   cluster sizes / RMSE';
for k = 1:numel(L)
    sz = arrayfun(@(c) c.n, L(k).CL);
    s1 = sprintf('%d ', sz);
    s2 = sprintf('%.3f ', L(k).rmse);
    lines{end+1} = sprintf('   T%-6d %6d   %8.4g   %8d   [%s] / [%s]', ...
        L(k).index, L(k).nData, L(k).radius, numel(L(k).CL), strtrim(s1), strtrim(s2)); %#ok<AGROW>
end
txt = strjoin(lines, sprintf('\n'));
if nargout == 0
    fprintf('%s\n', txt);
    clear txt
end
end
