function tf = so_runtime_exceeded(cfg)
%SO_RUNTIME_EXCEEDED Optional per-run wall-clock guard.
%
% Returns false unless cfg.runtime.startTime (a tic handle) and
% cfg.runtime.limitSeconds are both present.  Used by the ensemble so a
% single pathological source point cannot stall a 50-trial sweep.
tf = false;
if ~isfield(cfg, 'runtime'), return; end
if ~isfield(cfg.runtime, 'startTime') || ~isfield(cfg.runtime, 'limitSeconds'), return; end
if isempty(cfg.runtime.limitSeconds) || ~isfinite(cfg.runtime.limitSeconds), return; end
tf = toc(cfg.runtime.startTime) > cfg.runtime.limitSeconds;
end
