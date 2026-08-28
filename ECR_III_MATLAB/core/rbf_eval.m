function Y = rbf_eval(net, X)
%RBF_EVAL  Evaluate an RBF network trained by RBF_TRAIN.
%
%   Y = RBF_EVAL(NET,X) with X (d x M) returns Y (q x M).

M = size(X,2);
if isempty(net.W)
    Y = repmat(net.const, 1, M);
    return
end
Xn  = (X - repmat(net.mu,1,M)) ./ repmat(net.sd,1,M);
Phi = rbf_design(Xn, net.C, net.sig);
Y   = (Phi*net.W)';
end
