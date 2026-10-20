%===============================================================================
% Runs the recording protocols on the isolated channel of the modified Piezo1 model and measures
% each quantity the way the experiment did.
%
%   M = measure_Piezo1(X, D, all)
%
%   X   : options from ChR2_Piezo1_opts
%   D   : targets from Piezo1_targets
%   all : false computes the three fitted quantities only, true adds the checks
%
% PROTOCOLS (all at the recording temperature D.T_rec, from rest):
%   decay    whole cell poke: the stimulus rises in 1 ms to saturation and is
%            held for 300 ms at the test voltage. A single exponential with an
%            offset is fitted from the peak to the stimulus offset, as in Wu et
%            al. 2017.
%   deact    pressure clamp, -80 mV: 40 ms at saturation, then released within
%            2 ms. A single exponential with an offset is fitted to the 30 ms
%            after the release.
%   recovery pressure clamp, -80 mV: two steps of 100 ms, a delay between them.
%            Recovery = (peak 2 - end of step 1) / (peak 1 - end of step 1).
%   human    as decay, with a poke of 150 ms at -80 mV (Albuisson et al. 2013).
%   train    a 300 ms step, 10 s of rest, then 30 ms steps for 4 s at the given
%            rate. First and last peak as a fraction of the peak of the step.
%   sust     current at the end of the 300 ms step at -80 mV over its peak.
% Saturation is twice the reference tension of the activation gate.
%===============================================================================

function M = measure_Piezo1(X, D, all)

T    = D.T_rec;
Fsat = 2*X.F_ref;                 % kPa, opens the activation gate fully
rise = 0.001;                     % second, poke
ramp = 0.002;                     % second, pressure clamp

%% Fitted quantities
Vd          = unique(D.decay.V_mV,'stable');
tauV        = arrayfun(@(v) decay_tau(X,v,T,Fsat,rise,0.300),Vd);
M.decay_V   = Vd;
M.decay_tau = interp1(Vd,tauV,D.decay.V_mV);                 % ms, one value per data point
M.deact_tau = deact_tau(X,T,Fsat,ramp);                      % ms
M.half      = half_delay(X,T,Fsat,ramp);                     % second
if ~all
    return
end

%% Checks
M.human_tau = decay_tau(X,-80,T,Fsat,rise,0.150);            % ms
M.zheng80   = decay_tau(X,-80,T,Fsat,rise,0.300);            % ms
M.wufar_tau = arrayfun(@(v) decay_tau(X,v,T,Fsat,rise,0.300),D.wufar.V_mV);
M.rec       = arrayfun(@(d) recovery(X,T,Fsat,ramp,d),D.rec.delay);
[M.train_first,M.train_last] = arrayfun(@(f) train(X,T,Fsat,ramp,f),D.train.f_Hz);
S           = poke(X,-80,T,Fsat,rise,0.300);
M.sust      = S.g(end)/max(S.g);
M.peak_open = max(S.g);                                      % open fraction at the peak, -80 mV
end

%===============================================================================
% Local functions
%===============================================================================

function S = poke(X, V_mV, T, Fsat, rise, hold)
% One stimulus that rises in the given time and is held; output every 20 microseconds.
tb = [0 rise rise+hold];
S  = channel_run(X,V_mV,T,tb,[0 Fsat Fsat],(0:2e-5:rise+hold)');
end

function tau = decay_tau(X, V_mV, T, Fsat, rise, hold)
% Time constant in ms of the decay from the peak to the stimulus offset.
S      = poke(X,V_mV,T,Fsat,rise,hold);
[~,ip] = max(S.g);
tau    = 1e3*fit_exp(S.t(ip:end)-S.t(ip),S.g(ip:end));
end

function tau = deact_tau(X, T, Fsat, ramp)
% Time constant in ms of the decay after the release of a 40 ms step at -80 mV.
on  = 0.040;
tb  = [0 ramp on on+ramp on+ramp+0.030];
S   = channel_run(X,-80,T,tb,[0 Fsat Fsat 0 0],(on+ramp:2e-5:on+ramp+0.030)');
tau = 1e3*fit_exp(S.t-S.t(1),S.g);
end

function r = recovery(X, T, Fsat, ramp, delay)
% Two steps of 100 ms; delay from the end of the first to the start of the second.
d    = 0.100;
t2   = d + ramp + delay;
tb   = [0 ramp d d+ramp t2 t2+ramp t2+d];
Fb   = [0 Fsat Fsat 0 0 Fsat Fsat];
tout = [(0:2e-5:d)'; (t2:2e-5:t2+d)'];
S    = channel_run(X,-80,T,tb,Fb,tout);
w1   = S.t <= d;
p1   = max(S.g(w1));
e1   = S.g(find(w1,1,'last'));
p2   = max(S.g(~w1));
r    = (p2 - e1)/(p1 - e1);
end

function d = half_delay(X, T, Fsat, ramp)
% Delay in seconds at which the recovery is one half, by bisection on a logarithmic axis.
lo = -3;
hi = 1;
for k = 1:40
    mid = 0.5*(lo+hi);
    if recovery(X,T,Fsat,ramp,10^mid) < 0.5
        lo = mid;
    else
        hi = mid;
    end
end
d = 10^(0.5*(lo+hi));
end

function [first, last] = train(X, T, Fsat, ramp, f)
% A 300 ms step, 10 s of rest, then 30 ms steps for 4 s at f hertz.
d0   = 0.300;
w    = 0.030;
t0   = d0 + ramp + 10;
n    = round(4*f);
on   = t0 + (0:n-1)/f;
tb   = [0 ramp d0 d0+ramp];
Fb   = [0 Fsat Fsat 0];
for k = 1:n
    tb = [tb on(k) on(k)+ramp on(k)+w on(k)+w+ramp]; %#ok<AGROW>
    Fb = [Fb 0 Fsat Fsat 0];                          %#ok<AGROW>
end
tout = [(0:2e-5:d0)'; (on(1):2e-5:on(1)+w+ramp)'; (on(end):2e-5:on(end)+w+ramp)'];
tout = unique(tout);
S    = channel_run(X,-80,T,tb,Fb,tout);
pref = max(S.g(S.t <= d0));
first = max(S.g(S.t >= on(1)   & S.t <= on(1)+w+ramp))/pref;
last  = max(S.g(S.t >= on(end) & S.t <= on(end)+w+ramp))/pref;
end

function tau = fit_exp(t, y)
% Least squares fit of y = c + a*exp(-t/tau). For a given tau, c and a follow
% from a linear fit, so only tau is searched. Returns tau in the unit of t.
tau = fminbnd(@(x) resid(x,t,y),2e-4,2,optimset('TolX',1e-9));
end

function s = resid(tau, t, y)
A = [ones(size(t)) exp(-t/tau)];
c = A\y;
s = sum((y - A*c).^2);
end
