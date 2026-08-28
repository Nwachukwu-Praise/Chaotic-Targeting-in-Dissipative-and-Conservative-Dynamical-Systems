function [family, profile] = so_build_curve_family(component, direction, maxIter, cfg, profile, familyKind)
%SO_BUILD_CURVE_FAMILY Build all images/preimages once for diagonal reuse.
family = cell(1, maxIter + 1);
for n = 0:maxIter
    family{n + 1} = so_build_curve(component, direction, n, cfg);
end
if nargin >= 6
    switch familyKind
        case 'forward'
            profile.forwardFamilyBuilds = profile.forwardFamilyBuilds + 1;
        case 'backward'
            profile.backwardFamilyBuilds = profile.backwardFamilyBuilds + 1;
    end
end
profile = so_update_peak_curve_count(profile, family);
end

