function [cache, profile] = so_build_backward_cache(targetComponents, cfg, profile)
%SO_BUILD_BACKWARD_CACHE Build target preimage families once per stage.
cache.components = targetComponents;
cache.families = cell(1, numel(targetComponents));
for i = 1:numel(targetComponents)
    [cache.families{i}, profile] = so_build_curve_family(targetComponents{i}, -1, ...
        cfg.maxBackwardIterations, cfg, profile, 'backward');
end
end

