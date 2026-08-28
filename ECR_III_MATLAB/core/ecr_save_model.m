function ecr_save_model(fname, varargin)
%ECR_SAVE_MODEL  Save trained ECR models to a .mat file.
%
%   ECR_SAVE_MODEL(FNAME, M1, M2, ...) stores the models with the system's
%   function handles removed (S.step / S.obs / S.init / S.hfromz), because
%   Octave cannot serialise anonymous function handles.  Only S.name and the
%   numeric description survive; rebuild the live system with the matching
%   sys_*() constructor and put it back with
%
%       m = ECR_LOAD_MODEL('models.mat');
%       m{1}.S = sys_lorenz_real();

models = cell(1, numel(varargin));
for k = 1:numel(varargin)
    m = varargin{k};
    S = m.S;
    for f = {'step','obs','init','hfromz'}
        if isfield(S, f{1}), S.(f{1}) = []; end
    end
    if isfield(S,'cfg') && isstruct(S.cfg) && isfield(S.cfg,'meas')
        S.cfg.meas = [];
    end
    m.S = S;
    models{k} = m;
end
save(fname, 'models', '-v7');
end
