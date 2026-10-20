%===============================================================================
% Integrates hiMCES2026_ChR2_Piezo1 through light and mechanical pulses
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
%   o = run_segments(X, Y0, t0, t1, dt, classicVSoptic)
%
%   X              : options (ChR2_Piezo1_opts.m, OptoSono_settings.m); the light
%                    and mechanical pulse trains are defined in X
%   Y0             : initial state, 39 x 1 or 40 x 1
%   t0, t1         : start and end time, second
%   dt             : output interval, second
%   classicVSoptic : 0 standard solution, 1 Tyrode solution (host_arguments.m)
%
%   o              : t (s), Vm_mV, Cai_nM, Nai_mM, gates (states 34 and 35),
%                    tension_kPa (myofilament stress), I_Piezo1 (A/F), F_Piezo1
%                    (stress sensed by the Piezo1 activation gate, kPa), I_ChR2
%                    (A/F), Y (all states) and Yend (state at t1)
%
% Every pulse switches the right hand side discontinuously, so the integration
% is cut at every pulse edge. Each segment is integrated strictly inside its
% own regime, from a + 1e-9 to b - 1e-9 s, and the state is carried across the
% gap. No electrical stimulus is applied.
%===============================================================================

function o = run_segments(X, Y0, t0, t1, dt, classicVSoptic)

%% Pulse edges
edges = [t0 t1];
if X.optIrr > 0 && X.optN > 0
    on    = X.optT0 + (0:X.optN-1)/X.optFreq;
    edges = [edges on on+X.optDur];
end
if X.touchF > 0 && X.touchN > 0
    on    = X.touchT0 + (0:X.touchN-1)/X.touchFreq;
    edges = [edges on on+X.touchDur];
end
edges   = unique(round(edges(edges >= t0 & edges <= t1),9));

%% Integration, one segment at a time
A       = host_arguments(0,0.5,classicVSoptic);
options = odeset('MaxStep',5e-4,'InitialStep',1e-7,'RelTol',1e-7,'AbsTol',1e-10);
gap     = 1e-9;                          % second, distance kept from every edge
y       = Y0(:);
tout    = (t0:dt:t1)';
Y       = zeros(numel(tout),numel(y));
te      = tout;
for j = 1:numel(edges)-1
    a   = edges(j);
    b   = edges(j+1);
    sol = ode15s(@(t,x) hiMCES2026_ChR2_Piezo1(t,x,A{:},X),[a+gap b-gap],y,options);
    if sol.x(end) < b-gap-1e-12
        error('integration stopped at t = %.9g s',sol.x(end));
    end
    w = tout >= a & tout < b;
    if j == numel(edges)-1
        w = tout >= a & tout <= b;
    end
    if any(w)
        te(w)  = min(max(tout(w),a+gap),b-gap);
        Y(w,:) = deval(sol,te(w))';
    end
    y = sol.y(:,end);
end

%% Model outputs at the output times
D = zeros(numel(tout),4);
for i = 1:numel(tout)
    [~,q]  = hiMCES2026_ChR2_Piezo1(te(i),Y(i,:)',A{:},X);
    D(i,:) = q([21 49 51 56]);
end

o.t           = tout;
o.Vm_mV       = 1e3*Y(:,1);
o.Cai_nM      = 1e6*Y(:,3);
o.Nai_mM      = Y(:,18);
o.gates       = Y(:,[34 35]);
o.tension_kPa = D(:,1);
o.I_Piezo1    = D(:,2);              % A/F
o.F_Piezo1    = D(:,3);              % kPa, stress the activation gate senses
o.I_ChR2      = D(:,4);              % A/F
o.Y           = Y;
o.Yend        = y;
end
