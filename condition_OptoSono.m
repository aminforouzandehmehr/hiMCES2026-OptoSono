%===============================================================================
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
% PURPOSE:
%   Conditions the Fig. 8A cell under 0.5 Hz light pacing at 304 K in Tyrode.
%   The integration loop and convergence criterion are those used for the
%   manuscript's stored initial conditions.
%
%   st = condition_OptoSono('Piezo1')     modified Piezo1 model, 40 states
%   st = condition_OptoSono('Piezo1off')  Piezo1 current off, 39 states
%   With no argument the Piezo1 cell is conditioned.
%
% INITIAL CONDITIONS:
%   initial_conditions/start_31C_<cellType>.mat, variable Y
%
% OUTPUT:
%   results/settled_<cellType>_31C.mat, variable st: Y (final state), X
%   (options), T (conditioning time, s), ok (convergence flag), last (last
%   light-paced cycle). The supplied initial conditions are not overwritten.
%
% CONVERGENCE:
%   At most 20 blocks of 200 s. The relative change of APD90, peak calcium
%   and initial sodium between blocks must be below 2e-4, and the absolute
%   change of both Piezo1 gates below 1e-3. The first 99 cycles of each block
%   return endpoints; the final cycle is sampled every 2e-4 s.
%===============================================================================

function st = condition_OptoSono(cellType)

%% Housekeeping
if nargin == 0
    cellType = 'Piezo1';
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'model'));
addpath(fullfile(here,'analysis'));
outdir = fullfile(here,'results');
if ~isfolder(outdir)
    mkdir(outdir);
end

%% ============================================================================
%  1) PARAMETERS AND INITIAL CONDITIONS
% ============================================================================
S.X   = OptoSono_settings(cellType);
L     = load(fullfile(here,'initial_conditions',sprintf('start_31C_%s.mat',cellType)));
S.Y0  = L.Y(:);
S.cvo = 1;                               % Tyrode solution

%% ============================================================================
%  2) CONDITIONING
% ============================================================================
y    = S.Y0;
t    = 0;
prev = [];
gp   = [];
ok   = false;
S.X.optN = 1;
tic
for b = 1:20
    for c = 1:99
        S.X.optT0 = t;
        o = run_segments(S.X,y,t,t+2,2,S.cvo);
        y = o.Yend;
        t = t+2;
    end
    S.X.optT0 = t;
    o = run_segments(S.X,y,t,t+2,2e-4,S.cvo);
    y = o.Yend;
    t = t+2;
    u = upstroke_times(o.t,o.Vm_mV);
    if isempty(u)
        m = [nan max(o.Cai_nM) o.Nai_mM(1)];
    else
        m = [apd90_ms(o.t,o.Vm_mV,u(1)) max(o.Cai_nM) o.Nai_mM(1)];
    end
    g = o.gates(1,:);
    if ~isempty(prev)
        dr = max(abs(m./prev-1));
        dg = max(abs(g-gp));
        fprintf('%s, %g s: APD90 %.1f ms, Ca %.1f nM, Na_i %.3f mM, upstrokes %d, drift %.2g, gates %.2g\n', ...
                cellType,t,m,numel(u),dr,dg);
        if dr < 2e-4 && dg < 1e-3
            ok = true;
            break
        end
    end
    prev = m;
    gp   = g;
end

%% ============================================================================
%  3) OUTPUT
% ============================================================================
st.Y  = y;
st.X  = S.X;
st.T  = t;
st.ok = ok;
st.last = o;
save(fullfile(outdir,sprintf('settled_%s_31C.mat',cellType)),'st');
fprintf('%s settled %d after %g s (%.0f s)\n',cellType,ok,t,toc);
end
