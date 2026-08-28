function Phi = rbf_design(Xn, C, sig)
%RBF_DESIGN  Gaussian RBF design matrix with a bias column.
%
%   PHI = RBF_DESIGN(XN,C,SIG) with XN (d x n) standardised inputs, C (d x k)
%   centres and SIG (1 x k) widths returns PHI (n x (k+1)):
%
%       PHI(i,j) = exp(-||xn_i - c_j||^2 / (2*sig_j^2)),  PHI(:,k+1) = 1.

n = size(Xn,2);
k = size(C,2);
Phi = ones(n, k+1);
for j = 1:k
    d2 = sum((Xn - repmat(C(:,j),1,n)).^2, 1);
    Phi(:,j) = exp(-d2(:)/(2*sig(j)^2));
end
end
