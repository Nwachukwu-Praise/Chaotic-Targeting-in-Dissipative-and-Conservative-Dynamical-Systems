function tf = ecr_in_pi(P, S, tol)
%ECR_IN_PI  Test membership of the allowable parameter set Pi (Eq. 3).
%
%   tf = ECR_IN_PI(P,S) returns a 1 x M logical vector which is true where the
%   parameter vector P(:,k) satisfies
%
%       p_nom^i - dp_max^i < p^i < p_nom^i + dp_max^i ,   i = 1..r      (3)
%
%   tf = ECR_IN_PI(P,S,TOL) allows a relative tolerance (default 1e-9).

if nargin < 3 || isempty(tol), tol = 1e-9; end
D  = abs(P - repmat(S.pnom, 1, size(P,2)));
tf = all(D <= repmat(S.dpmax*(1+tol), 1, size(P,2)), 1);
tf = tf & all(isfinite(P), 1);
end
