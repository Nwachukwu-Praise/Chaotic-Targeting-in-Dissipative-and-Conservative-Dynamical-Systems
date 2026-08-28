function [isValid, differences] = validate_bisection_noise_config( ...
    storedConfig, currentConfig)
%VALIDATE_BISECTION_NOISE_CONFIG Field-by-field configuration comparison.

isValid = true;
differences = strings(0, 1);
if ~isstruct(storedConfig) || ~isstruct(currentConfig)
    isValid = false;
    differences(end + 1, 1) = "configuration is not a struct";
    return;
end

fields = fieldnames(currentConfig);
for i = 1:numel(fields)
    fieldName = fields{i};
    if ~isfield(storedConfig, fieldName)
        isValid = false;
        differences(end + 1, 1) = "missing field: " + fieldName;
        continue;
    end
    if ~values_match(storedConfig.(fieldName), currentConfig.(fieldName))
        isValid = false;
        differences(end + 1, 1) = "different field: " + fieldName;
    end
end
end

function tf = values_match(a, b)
if isnumeric(a) || islogical(a)
    tf = (isnumeric(b) || islogical(b)) && isequaln(size(a), size(b)) && ...
        all(abs(a(:) - b(:)) <= 100 * eps(max(1, max(abs([a(:); b(:)]))))) ;
elseif ischar(a)
    tf = ischar(b) && strcmp(a, b);
elseif isstring(a)
    tf = isstring(b) && isequal(a, b);
else
    tf = isequaln(a, b);
end
end
