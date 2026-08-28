function Dm = ecr_nmd(CL, Z)
%ECR_NMD  Normalised Mahalanobis Distance of states to clusters (Eq. 8).
%
%   DM = ECR_NMD(CL,Z) returns DM (nc x M) with
%
%       DM(j,k) = sqrt( (z_k - mu_j)' * N_j * (z_k - mu_j) )               (8)
%
%   i.e. the NMD (not its square) of every state to every cluster.  CL is the
%   cluster array built by ECR_CLUSTER_STATS.

nc = numel(CL);
M  = size(Z,2);
Dm = zeros(nc, M);
for j = 1:nc
    Xc = Z - repmat(CL(j).mu, 1, M);
    Dm(j,:) = sqrt(max(sum(Xc .* (CL(j).Ninv*Xc), 1), 0));
end
end
