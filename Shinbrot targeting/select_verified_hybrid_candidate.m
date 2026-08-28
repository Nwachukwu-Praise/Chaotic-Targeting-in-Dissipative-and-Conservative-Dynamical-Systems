function [selected, uniqueCandidates] = ...
    select_verified_hybrid_candidate(candidates)
%SELECT_VERIFIED_HYBRID_CANDIDATE Deduplicate, then rank by |p| and error.

if isempty(candidates)
    error('At least one verified candidate is required.');
end

pValues = [candidates.p].';
errors = [candidates.error].';
if any(~isfinite(pValues)) || any(~isfinite(errors))
    error('Verified candidates must have finite p and error values.');
end

dedupTolerance = 1e-12;
[~, order] = sortrows([abs(pValues), errors, pValues], [1, 2, 3]);
ordered = candidates(order);
keep = true(size(ordered));
for i = 2:numel(ordered)
    previousP = [ordered(keep(1:i - 1)).p];
    keep(i) = all(abs(ordered(i).p - previousP) > dedupTolerance);
end
uniqueCandidates = ordered(keep);
selected = uniqueCandidates(1);
end
