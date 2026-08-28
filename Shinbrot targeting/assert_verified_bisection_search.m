function assert_verified_bisection_search(search, verifiedIdentifier, context, target)
%ASSERT_VERIFIED_BISECTION_SEARCH Validate Shinbrot bisection provenance.

if nargin < 2 || isempty(verifiedIdentifier)
    verifiedIdentifier = verified_bisection_identifier();
end
if nargin < 3 || isempty(context)
    context = 'search result';
end
if nargin < 4
    target = struct();
end

if ~isstruct(search)
    error('%s is not a search result struct.', context);
end
if ~isfield(search, 'method') || ~strcmp(search.method, verifiedIdentifier)
    error('%s did not use dispatcher %s.', context, verifiedIdentifier);
end
if strcmp(verifiedIdentifier, 'shinbrotPaperBisection')
    if ~isfield(search, 'paperAlgorithmPath') || ~search.paperAlgorithmPath
        error('%s is not marked as the paper-bisection algorithm path.', context);
    end
end
if isfield(search, 'gridEvaluations') && search.gridEvaluations ~= 0
    error('%s contains hybrid-style full-grid evaluations.', context);
end
if isfield(search, 'crossingOrderTestUsed') && ~search.crossingOrderTestUsed
    error('%s did not record the discontinuity-aware crossing-order test.', context);
end
if isfield(search, 'selectionMethod') && contains(search.selectionMethod, 'grid-first')
    error('%s used hybrid grid-first candidate selection.', context);
end
if isfield(search, 'directGridHitsByHorizon') && any(search.directGridHitsByHorizon ~= 0)
    error('%s contains hybrid direct-grid-hit metadata.', context);
end
if isfield(search, 'adjacentBracketsByHorizon') && any(search.adjacentBracketsByHorizon ~= 0)
    error('%s contains hybrid adjacent-grid-bracket metadata.', context);
end
if strcmp(verifiedIdentifier, 'shinbrotPaperBisection') && ...
        isfield(search, 'found') && search.found
    if ~isfield(search, 'paperReplayVerified') || ~search.paperReplayVerified
        error('%s lacks the paper direct replay consistency check.', context);
    end
    if isfield(target, 'tolerance') && isfield(search, 'finalTargetError') && ...
            ~(isfinite(search.finalTargetError) && ...
            search.finalTargetError <= target.tolerance)
        error('%s replay target error exceeds the declared tolerance.', context);
    end
    if isfield(search, 'refinementsCompleted') && ...
            isfield(search, 'finalBracketBisectionIterations') && ...
            search.finalBracketBisectionIterations > 0 && ...
            search.refinementsCompleted ~= 24
        error('%s did not complete 24 paper final-bracket refinements.', context);
    end
end
end
