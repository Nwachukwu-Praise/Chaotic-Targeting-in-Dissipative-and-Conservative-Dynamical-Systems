function D = ecr_cache_data(S, opt, fname)
%ECR_CACHE_DATA  Generate the training triples once and reuse them.
%
%   D = ECR_CACHE_DATA(S,OPT,FNAME) loads FNAME if it exists and its stored
%   options match OPT (same seed / size), otherwise it runs
%   ECR_GENERATE_DATA and saves the result.  Data generation is by far the
%   most expensive part of a Lorenz experiment, so caching keeps the demos
%   quick to re-run.

if nargin < 3 || isempty(fname)
    fname = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                     'results', sprintf('data_%s.mat', S.name));
end
d = fileparts(fname);
if ~isempty(d) && exist(d, 'dir') ~= 7, mkdir(d); end

if exist(fname, 'file') == 2
    L = load(fname);
    if isfield(L,'D') && isfield(L.D,'opt') && ...
       L.D.opt.nData == opt.nData && L.D.opt.seed == opt.seed && ...
       L.D.opt.trajLen == opt.trajLen && size(L.D.Z,1) == S.N
        D = L.D;
        if opt.verbose
            fprintf('[data] reusing %s (%d triples)\n', fname, size(D.Z,2));
        end
        return
    end
end

D = ecr_generate_data(S, opt);
save(fname, 'D', '-v7');
if opt.verbose, fprintf('[data] saved %s\n', fname); end
end
