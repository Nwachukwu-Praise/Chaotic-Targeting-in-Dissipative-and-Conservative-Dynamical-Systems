function opts = so_parse_options(args, defaults)
%SO_PARSE_OPTIONS Minimal name/value parser (no toolbox dependency).
%
% inputParser lives in base MATLAB but is slow and verbose for this; this
% helper keeps the runners readable and works unchanged in Octave.
opts = defaults;
if isempty(args)
    return;
end
if mod(numel(args), 2) ~= 0
    error('SchroerOtt:BadOptions', 'Options must be name/value pairs.');
end
known = fieldnames(defaults);
for i = 1:2:numel(args)
    name = args{i};
    if ~(ischar(name) || isstring(name))
        error('SchroerOtt:BadOptionName', 'Option names must be text.');
    end
    name = char(name);
    match = known(strcmpi(known, name));
    if isempty(match)
        error('SchroerOtt:UnknownOption', 'Unknown option "%s". Known: %s.', ...
            name, strjoin(known.', ', '));
    end
    opts.(match{1}) = args{i + 1};
end
end
