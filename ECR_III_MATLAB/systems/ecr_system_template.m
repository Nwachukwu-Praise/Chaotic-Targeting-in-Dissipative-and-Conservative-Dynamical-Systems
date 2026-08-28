function S = ecr_system_template()
%ECR_SYSTEM_TEMPLATE  Empty system description used by every sys_*.m file.
%
%   Every chaotic system handled by this ECR toolbox is described by a struct
%   with the fields below.  The struct hides the difference between
%
%     (a) plain discrete maps            z_{n+1} = G(z_n,p_n)          (Eq. 1)
%     (b) Poincare maps of flows         (state on the surface of section)
%     (c) delay-coordinate observations  z_n = [s(t_n) s(t_n-T) ...]
%
%   by separating the *hidden* state H (everything needed to propagate the
%   system) from the *observed* state Z (what the controller is allowed to
%   see, i.e. the vector z_n of Eq. 1).  For plain maps H == Z.
%
%   Fields
%     name     char, short identifier
%     N        dimension of the observed state z            (N x 1)
%     nh       dimension of the hidden state h
%     r        number of control parameters                  (r x 1)
%     zstar    target (unstable fixed point of the map)       N x 1   (Eq. 2)
%     pnom     nominal parameter vector                       r x 1
%     dpmax    maximum allowable parameter change dp_max      r x 1   (Eq. 3)
%     delta    radius of the OGY / S0 region                  scalar  (Eq. 4)
%     rms      RMS of the measured state (used to scale noise)
%     pnames   1 x r cell array of parameter names
%     type     'map' | 'poincare' | 'delay'
%     step     @(H,P) -> [Hnext, Znext]   propagate one step (columns = samples)
%     obs      @(H)   -> Z                observed state of a hidden state
%     init     @(m)   -> H                m random hidden states on the attractor
%     hfromz   @(Z)   -> H  or []         lift an observed state back to a hidden
%                                         state (only exists for 'map' systems)
%     notes    char, free text
%
%   See also SYS_LOGISTIC, SYS_HENON, SYS_LORENZ_REAL, SYS_LORENZ_DELAY,
%            SYS_DUFFING_DELAY.

S = struct( ...
    'name',   '', ...
    'N',      0,  ...
    'nh',     0,  ...
    'r',      0,  ...
    'zstar',  [], ...
    'pnom',   [], ...
    'dpmax',  [], ...
    'delta',  0,  ...
    'rms',    1,  ...
    'pnames', {{}}, ...
    'type',   'map', ...
    'step',   [], ...
    'obs',    [], ...
    'init',   [], ...
    'hfromz', [], ...
    'notes',  '');
end
