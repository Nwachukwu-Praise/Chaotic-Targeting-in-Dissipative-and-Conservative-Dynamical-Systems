function v = so_interquartile_range(x)
%SO_INTERQUARTILE_RANGE IQR without the Statistics Toolbox.
%
% iqr and prctile live in the Statistics and Machine Learning Toolbox, which
% this project does not otherwise need.  Linear-interpolated quartiles on
% sample points at (i-0.5)/n reproduce the default prctile convention.
x = sort(x(:));
n = numel(x);
if n == 0
    v = NaN;
elseif n == 1
    v = 0;
else
    v = percentile(x, 75) - percentile(x, 25);
end
end

function q = percentile(sortedX, pct)
n = numel(sortedX);
pos = pct / 100 * n + 0.5;
if pos <= 1
    q = sortedX(1);
elseif pos >= n
    q = sortedX(n);
else
    lo = floor(pos);
    frac = pos - lo;
    q = (1 - frac) * sortedX(lo) + frac * sortedX(lo + 1);
end
end
