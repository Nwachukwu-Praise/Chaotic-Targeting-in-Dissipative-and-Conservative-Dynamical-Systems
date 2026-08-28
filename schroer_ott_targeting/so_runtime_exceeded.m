function tf = so_runtime_exceeded(cfg)
%SO_RUNTIME_EXCEEDED Optional wall-clock guard for long comparison runs.
tf = false;
if ~isfield(cfg, 'runtime') || ~isfield(cfg.runtime, 'limitSeconds') || ...
        ~isfield(cfg.runtime, 'startTime')
    return;
end
if isempty(cfg.runtime.limitSeconds) || ~isfinite(cfg.runtime.limitSeconds) || cfg.runtime.limitSeconds <= 0
    return;
end
tf = toc(cfg.runtime.startTime) > cfg.runtime.limitSeconds;
end
