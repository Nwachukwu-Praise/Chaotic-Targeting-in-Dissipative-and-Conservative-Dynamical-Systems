function classification = so_classify_trace(traceValue, traceTolerance)
%SO_CLASSIFY_TRACE Hyperbolic/elliptic/parabolic trace classification.
if traceValue > 2 + traceTolerance
    classification = 'direct-hyperbolic';
elseif traceValue < -2 - traceTolerance
    classification = 'inverse-hyperbolic';
elseif abs(traceValue) < 2 - traceTolerance
    classification = 'elliptic';
else
    classification = 'near-parabolic';
end
end

