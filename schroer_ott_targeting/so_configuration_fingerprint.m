function fingerprint = so_configuration_fingerprint(value)
%SO_CONFIGURATION_FINGERPRINT Stable lightweight fingerprint for settings.
canonical = canonical_text(value);
bytes = uint8(char(canonical));
hash = uint32(2166136261);
prime = uint32(16777619);
for i = 1:numel(bytes)
    hash = uint32(bitxor(hash, uint32(bytes(i))));
    hash = uint32(mod(uint64(hash) * uint64(prime), uint64(2)^32));
end
fingerprint = string(sprintf('fnv1a32:%08x', hash));
end

function txt = canonical_text(value)
if isstruct(value)
    if numel(value) ~= 1
        parts = strings(1, numel(value));
        for i = 1:numel(value)
            parts(i) = canonical_text(value(i));
        end
        txt = "[" + strjoin(parts, ",") + "]";
        return;
    end
    names = sort(string(fieldnames(value)));
    parts = strings(1, numel(names));
    for i = 1:numel(names)
        name = char(names(i));
        parts(i) = names(i) + ":" + canonical_text(value.(name));
    end
    txt = "{" + strjoin(parts, ",") + "}";
elseif istable(value)
    txt = evalc('disp(value)');
elseif iscell(value)
    parts = strings(1, numel(value));
    for i = 1:numel(value)
        parts(i) = canonical_text(value{i});
    end
    txt = "{" + strjoin(parts, ",") + "}";
elseif isnumeric(value)
    txt = mat2str(value, 17);
elseif islogical(value)
    txt = mat2str(value);
elseif isstring(value)
    txt = "[""" + strjoin(value(:).', """,""") + """]";
elseif ischar(value)
    txt = string(value);
else
    txt = string(evalc('disp(value)'));
end
txt = string(txt);
end
