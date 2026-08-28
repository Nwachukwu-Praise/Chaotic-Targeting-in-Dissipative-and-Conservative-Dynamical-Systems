function tf = so_should_prune_probe(j, bestObjective)
%SO_SHOULD_PRUNE_PROBE Exact switch-probe pruning rule.
tf = isfinite(bestObjective) && j >= bestObjective;
end

