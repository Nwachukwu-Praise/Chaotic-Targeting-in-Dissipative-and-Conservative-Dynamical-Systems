function info = shinbrot_dependency_fingerprint(kind)
%SHINBROT_DEPENDENCY_FINGERPRINT Hash source files used by a workflow.
%
% kind is 'deterministic', 'noise' or 'all'.  The hashes make saved result
% validation sensitive to implementation edits, not just configuration text.

if nargin < 1 || isempty(kind)
    kind = 'all';
end
kind = char(kind);

deterministicFiles = { ...
    'search_parameter_to_target', ...
    'run_shinbrot_paper_bisection', ...
    'evaluate_shinbrot_crossing_signature', ...
    'next_valid_section_crossing', ...
    'poincare_event', ...
    'lorenz_rhs'};

noiseFiles = { ...
    'run_noisy_targeting_trial', ...
    'run_noisy_trial_paper_cycle_cadence', ...
    'run_noisy_trial_paper40_step', ...
    'run_noisy_trial_section_retarget', ...
    'run_noise_sweep', ...
    'select_noise_parallel_workers', ...
    'rk4_lorenz_step', ...
    'next_valid_section_crossing_noisy_rk4', ...
    'integrate_noisy_lorenz_fixed_step'};

switch lower(kind)
    case 'deterministic'
        functionNames = deterministicFiles;
    case 'noise'
        functionNames = noiseFiles;
    case 'all'
        functionNames = [deterministicFiles, noiseFiles];
    otherwise
        error('Unknown fingerprint kind: %s', kind);
end

componentNames = functionNames(:);
files = cell(size(componentNames));
sha256 = cell(size(componentNames));
componentText = cell(size(componentNames));

for i = 1:numel(componentNames)
    fileName = which(componentNames{i});
    if isempty(fileName)
        error('Required function is not on the MATLAB path: %s', ...
            componentNames{i});
    end
    files{i} = fileName;
    sourceText = normalize_newlines(fileread(fileName));
    sha256{i} = local_digest(sourceText, 'SHA-256');
    componentText{i} = sprintf('%s|%s|%s', ...
        componentNames{i}, fileName, sha256{i});
end

info.kind = kind;
info.componentNames = componentNames;
info.files = files;
info.sha256 = sha256;
info.signatureText = strjoin(componentText, newline);
info.aggregateSHA256 = local_digest(info.signatureText, 'SHA-256');
end

function text = normalize_newlines(text)
text = strrep(text, sprintf('\r\n'), newline);
text = strrep(text, sprintf('\r'), newline);
end

function digest = local_digest(text, algorithm)
bytes = uint8(text);
try
    engine = java.security.MessageDigest.getInstance(algorithm);
    engine.update(typecast(bytes, 'int8'));
    raw = typecast(engine.digest(), 'uint8');
    digest = lower(reshape(dec2hex(raw, 2).', 1, []));
catch
    weights = mod(1:numel(bytes), 251) + 1;
    digest = sprintf('nojvm%012.0f', ...
        mod(sum(double(bytes) .* weights), 1e12));
end
end
