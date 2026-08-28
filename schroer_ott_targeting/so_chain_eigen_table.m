function tbl = so_chain_eigen_table(chains, includeElliptic)
%SO_CHAIN_EIGEN_TABLE Monodromy spectrum of each retained periodic chain.
%
%   tbl = SO_CHAIN_EIGEN_TABLE(chains) returns one row per chain with the
%   trace-based classification and the stable/unstable eigenvalues of the
%   monodromy matrix DF^p(z0).
%
%   For an area-preserving map det DF^p = 1, so the two eigenvalues are
%   reciprocal.  The trace therefore classifies the chain completely:
%
%       trace >  2   direct-hyperbolic   (real, both positive)
%       trace < -2   inverse-hyperbolic  (real, both negative)
%      |trace| <  2  elliptic            (complex conjugate, unit modulus)
%       near +-2     near-parabolic
%
%   By default elliptic chains are omitted, because their eigenvalues are
%   complex and no stable/unstable splitting exists.  Pass includeElliptic
%   as true to list them with their complex pair.

if nargin < 2 || isempty(includeElliptic)
    includeElliptic = false;
end

if ~includeElliptic
    keep = ~strcmp({chains.classification}, 'elliptic');
    chains = chains(keep);
end

if isempty(chains)
    tbl = table('Size', [0 11], ...
        'VariableTypes', {'string','double','double','double','string','double', ...
        'double','double','double','double','double'}, ...
        'VariableNames', {'ChainID','Period','Winding','RotationNumber','Classification', ...
        'Trace','StableEigenvalue','UnstableEigenvalue','AbsStableEigenvalue', ...
        'AbsUnstableEigenvalue','LyapunovExponentPerIterate'});
    return;
end

n = numel(chains);
chainID = strings(n, 1);
period = zeros(n, 1);
winding = zeros(n, 1);
rotationNumber = zeros(n, 1);
classification = strings(n, 1);
traceValue = zeros(n, 1);
stableEigenvalue = complex(zeros(n, 1));
unstableEigenvalue = complex(zeros(n, 1));
absStable = zeros(n, 1);
absUnstable = zeros(n, 1);
lyapunov = zeros(n, 1);

for i = 1:n
    c = chains(i);
    e = local_eigen(c);
    chainID(i) = string(c.id);
    period(i) = c.period;
    winding(i) = c.winding;
    rotationNumber(i) = c.omega;
    classification(i) = string(c.classification);
    traceValue(i) = c.trace;
    stableEigenvalue(i) = e.stableEigenvalue;
    unstableEigenvalue(i) = e.unstableEigenvalue;
    absStable(i) = e.absStableEigenvalue;
    absUnstable(i) = e.absUnstableEigenvalue;
    lyapunov(i) = e.lyapunovExponentPerIterate;
end

% Present real spectra as real numbers; only elliptic rows stay complex.
if all(imag(stableEigenvalue) == 0) && all(imag(unstableEigenvalue) == 0)
    stableEigenvalue = real(stableEigenvalue);
    unstableEigenvalue = real(unstableEigenvalue);
end

tbl = table(chainID, period, winding, rotationNumber, classification, traceValue, ...
    stableEigenvalue, unstableEigenvalue, absStable, absUnstable, lyapunov, ...
    'VariableNames', {'ChainID','Period','Winding','RotationNumber','Classification', ...
    'Trace','StableEigenvalue','UnstableEigenvalue','AbsStableEigenvalue', ...
    'AbsUnstableEigenvalue','LyapunovExponentPerIterate'});
tbl = sortrows(tbl, {'RotationNumber','Period','ChainID'});
end

function e = local_eigen(chain)
%LOCAL_EIGEN Read cached eigen data, or recover it from the trace.
%
%   Chains produced by the current so_periodic_chain_from_root carry an
%   eigen substructure.  Chains loaded from an older stored catalogue do
%   not, so the reciprocal pair is recomputed here from the trace alone.
if isfield(chain, 'eigen') && ~isempty(chain.eigen)
    e = chain.eigen;
    return;
end
tr = chain.trace;
disc = tr^2 - 4;
if disc > 0
    root = sqrt(disc);
    a = 0.5 * (tr + root);
    b = 0.5 * (tr - root);
    if abs(a) >= abs(b)
        e.unstableEigenvalue = a;
        e.stableEigenvalue = b;
    else
        e.unstableEigenvalue = b;
        e.stableEigenvalue = a;
    end
else
    a = 0.5 * (tr + 1i * sqrt(max(-disc, 0)));
    e.unstableEigenvalue = a;
    e.stableEigenvalue = conj(a);
end
e.absStableEigenvalue = abs(e.stableEigenvalue);
e.absUnstableEigenvalue = abs(e.unstableEigenvalue);
e.lyapunovExponentPerIterate = log(e.absUnstableEigenvalue) / chain.period;
end
