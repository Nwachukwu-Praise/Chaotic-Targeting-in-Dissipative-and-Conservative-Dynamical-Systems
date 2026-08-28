function [fingerprint, info] = verified_bisection_implementation()
%VERIFIED_BISECTION_IMPLEMENTATION Fingerprint of the dispatched search.
%
%   VERIFIED_BISECTION_IDENTIFIER names the dispatcher case.  It does not
%   say which function that case calls, and it does not change when the
%   function behind it changes.  Saved noise and ensemble results were
%   therefore validated against configuration settings alone, which cannot
%   detect a change of implementation, and the stored results silently
%   remained valid while the search that produced them was replaced.  That
%   is how Tables 2 to 4 came to describe three different searches.
%
%   This function resolves the implementation actually dispatched by reading
%   search_parameter_to_target, and returns a fingerprint combining the
%   function name, its file size and an MD5 digest of its source.  Any edit
%   to the search changes the fingerprint, so any stored result validated
%   against it is correctly rejected rather than silently reused.
%
%   fingerprint : char, of the form
%                 'run_shinbrot_coverage_bisection|bytes=28114|md5=1a2b...'
%   info        : struct with the separate fields.
%
%   The implementation is not hard-coded here.  Repointing the dispatcher is
%   sufficient; this function follows it.

identifier = verified_bisection_identifier();

dispatcherFile = which('search_parameter_to_target');
if isempty(dispatcherFile)
    error('verified_bisection_implementation:DispatcherNotFound', ...
        'search_parameter_to_target is not on the MATLAB path.');
end
dispatcherText = fileread(dispatcherFile);

caseToken = sprintf('case ''%s''', identifier);
caseStart = strfind(dispatcherText, caseToken);
if isempty(caseStart)
    error('verified_bisection_implementation:CaseNotFound', ...
        'No dispatcher case for %s was found in %s.', identifier, dispatcherFile);
end
caseStart = caseStart(1) + numel(caseToken);

% The case block ends at the next case label or at otherwise.
remainder = dispatcherText(caseStart:end);
stopTokens = [strfind(remainder, 'case '), strfind(remainder, 'otherwise')];
if isempty(stopTokens)
    blockText = remainder;
else
    blockText = remainder(1:min(stopTokens) - 1);
end

% Ignore commentary; the call is the first non-comment run_* token.
lines = strsplit(blockText, newline);
codeLines = lines(~startsWith(strtrim(lines), '%'));
blockCode = strjoin(codeLines, newline);

called = regexp(blockCode, '(run_\w+)\s*\(', 'tokens', 'once');
if isempty(called)
    error('verified_bisection_implementation:CallNotFound', ...
        'No run_* call was found in the %s dispatcher case.', identifier);
end
functionName = called{1};

implementationFile = which(functionName);
if isempty(implementationFile)
    error('verified_bisection_implementation:ImplementationNotFound', ...
        'The dispatcher calls %s, which is not on the MATLAB path.', functionName);
end

fileInfo = dir(implementationFile);
dependencyInfo = shinbrot_dependency_fingerprint('deterministic');

info.identifier = identifier;
info.functionName = functionName;
info.file = implementationFile;
info.bytes = fileInfo.bytes;
info.dependencyFingerprint = dependencyInfo;
info.digest = dependencyInfo.aggregateSHA256;

fingerprint = sprintf('%s|bytes=%d|deterministic-sha256=%s', ...
    info.functionName, info.bytes, info.digest);
info.fingerprint = fingerprint;
end

function digest = local_digest(text)
%LOCAL_DIGEST MD5 of the source text, with a portable fallback.
%
%   Line endings are normalised first so that the fingerprint does not
%   change merely because a file was checked out on a different platform.
text = strrep(text, sprintf('\r\n'), newline);
text = strrep(text, sprintf('\r'), newline);
bytes = uint8(text);

try
    engine = java.security.MessageDigest.getInstance('MD5');
    engine.update(typecast(bytes, 'int8'));
    raw = typecast(engine.digest(), 'uint8');
    digest = lower(reshape(dec2hex(raw, 2).', 1, []));
catch
    % No JVM.  A weaker but still edit-sensitive checksum.
    weights = mod(1:numel(bytes), 251) + 1;
    digest = sprintf('nojvm%012.0f', mod(sum(double(bytes) .* weights), 1e12));
end
end
