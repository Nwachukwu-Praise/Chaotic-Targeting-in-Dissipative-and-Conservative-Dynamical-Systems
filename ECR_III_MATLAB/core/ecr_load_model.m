function models = ecr_load_model(fname, S)
%ECR_LOAD_MODEL  Load models written by ECR_SAVE_MODEL.
%
%   MODELS = ECR_LOAD_MODEL(FNAME)      cell array of models (no live system)
%   MODELS = ECR_LOAD_MODEL(FNAME,S)    re-attaches the live system S so that
%                                       ECR_CONTROL can be used again.

L = load(fname);
models = L.models;
if nargin >= 2 && ~isempty(S)
    for k = 1:numel(models)
        models{k}.S = S;
    end
end
end
