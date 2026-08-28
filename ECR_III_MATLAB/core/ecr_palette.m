function C = ecr_palette(n)
%ECR_PALETTE  Fixed categorical colours used by every figure of the toolbox.
%
%   C = ECR_PALETTE(N) returns an N x 3 RGB matrix.  The hues are assigned in
%   a fixed order (never cycled through a colormap), are separated for the
%   common colour-vision deficiencies, and all reach at least 3:1 contrast
%   against a white figure background, so a reader can tell control region
%   T_1 from T_2 in print, on screen and in greyscale-plus-legend.

base = [ 43 107 214; ...   % blue
        217  95   2; ...   % orange
         27 158 119; ...   % green
        117 112 179; ...   % purple
        140 109   0; ...   % dark gold
        176  48  96] / 255;% magenta
if nargin < 1 || isempty(n), n = size(base,1); end
C = base(mod(0:n-1, size(base,1)) + 1, :);
end
