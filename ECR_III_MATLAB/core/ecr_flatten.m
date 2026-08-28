function F = ecr_flatten(levels)
%ECR_FLATTEN  Pool every cluster C_ij of every control region T_i.
%
%   F = ECR_FLATTEN(LEVELS) returns
%     F.CL     1 x nc struct array of all clusters (ECR_CLUSTER_STATS format)
%     F.nets   1 x nc cell array of the matching networks NN_ij
%     F.level  1 x nc region index i of every cluster
%     F.jj     1 x nc cluster index j inside its region
%
%   ECR-III selects the region of the current state as the cluster with the
%   smallest Normalised Mahalanobis Distance over this pooled list
%   (Section 3.3).

CLs = {};  nets = {};  lev = [];  jj = [];
for k = 1:numel(levels)
    L = levels(k);
    for j = 1:numel(L.CL)
        CLs{end+1}  = L.CL(j);      %#ok<AGROW>
        nets{end+1} = L.nets{j};    %#ok<AGROW>
        lev(end+1)  = L.index;      %#ok<AGROW>
        jj(end+1)   = j;            %#ok<AGROW>
    end
end
if isempty(CLs)
    F = struct('CL', [], 'nets', {{}}, 'level', [], 'jj', []);
    return
end
CL = CLs{1};
for k = 2:numel(CLs), CL(k) = CLs{k}; end
F = struct('CL', CL, 'nets', {nets}, 'level', lev, 'jj', jj);
end
