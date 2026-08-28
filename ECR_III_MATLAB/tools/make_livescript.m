function out = make_livescript(src, dst)
%MAKE_LIVESCRIPT  Convert a sectioned .m script into a MATLAB live script.
%
%   MAKE_LIVESCRIPT                       converts demos/ECR_III_walkthrough.m
%                                         into demos/ECR_III_walkthrough.mlx
%   MAKE_LIVESCRIPT(SRC)                  converts SRC next to itself
%   MAKE_LIVESCRIPT(SRC,DST)              converts SRC into DST
%   OUT = MAKE_LIVESCRIPT(...)            returns the path that was written
%
%   Why this exists.  A .mlx file is a binary package that only MATLAB can
%   write, so the walkthrough ships as a .m file carrying live-script markup
%   (section titles, formatted text, LaTeX equations).  This function asks the
%   Live Editor to do the conversion, which takes a second and produces a
%   fully formatted live script with headings and rendered equations.
%
%   Requires MATLAB R2016a or newer.  If the internal Live Editor API is not
%   available in your release, open the .m file in MATLAB and use
%   Save > Save As... > MATLAB Live Code File (*.mlx) - the result is
%   identical.
%
%   Example
%       >> startup_ecr
%       >> make_livescript
%       >> open ECR_III_walkthrough.mlx

here = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1 || isempty(src)
    src = fullfile(here, 'demos', 'ECR_III_walkthrough.m');
end
if nargin < 2 || isempty(dst)
    [p, n] = fileparts(src);
    dst = fullfile(p, [n '.mlx']);
end

if exist(src, 'file') ~= 2
    error('make_livescript:noSource', 'cannot find %s', src);
end
if exist('OCTAVE_VERSION', 'builtin')
    error('make_livescript:octave', ...
        ['Live scripts are a MATLAB feature; Octave cannot write .mlx.  ' ...
         'Run %s directly instead - it is an ordinary script.'], src);
end

fprintf('converting\n  %s\n->\n  %s\n', src, dst);
try
    matlab.internal.liveeditor.openAndSave(char(src), char(dst));
catch err
    error('make_livescript:failed', ...
        ['The Live Editor conversion API is not available in this release ' ...
         '(%s).\nOpen %s in MATLAB and use Save > Save As... > MATLAB Live ' ...
         'Code File (*.mlx).'], err.message, src);
end

if exist(dst, 'file') ~= 2
    error('make_livescript:noOutput', 'conversion produced no file');
end
fprintf('done: %s\n', dst);
if nargout == 0
    clear out
else
    out = dst;
end
end
