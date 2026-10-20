%===============================================================================
% Fits the three gate parameters of the modified Piezo1 model to measured channel responses.
%
%   fit_Piezo1()
%
% WHAT IS FITTED:
%   tau_act   activation and deactivation time constant of the activation gate
%   k_inact   scale on the published voltage function of the gate time constant,
%             used here as the inactivation time constant of an open channel
%   tau_rec   recovery from inactivation
%   to twelve measured values (Piezo1_targets.m): the decay time constant of the
%   current at five voltages in two data sets, the deactivation time constant,
%   and the delay of half recovery. Each quantity is measured on the model with
%   the protocol of its experiment (measure_Piezo1.m). The cell model is not used.
%
% WHAT IS NOT FITTED:
%   The shape of the voltage function, the activation curve, the conductance and
%   the reversal potential stay as published. No cell level quantity (action
%   potential, calcium, tension, capture) enters the fit.
%
% RULE, written before the fit:
%   A fitted decay value passes if it lies within the plotted error bar of its
%   data point. Deactivation and half recovery pass within 20 %. A check passes
%   within two standard errors (or two standard deviations where that is what
%   the source gives), or within 20 % where the source gives no error.
%   Every target and every check stays in the table, passed or not.
%
% EXPECTATION, written before the fit:
%   Fitted: all twelve pass, with k_inact between 0.6 and 1.0 (the fitted decay
%   is slower than the gate time constant because activation is not instant).
%   Checks: human inactivation at -80 mV (8.6 ms) passes or is close. Zheng's
%   11.9 ms at -80 mV is missed (the function gives about 8 to 9 ms there).
%   Wu's values at +60 mV and above are missed, the model being too fast.
%   The recovery curve is missed at 10 to 40 ms (model too slow) and from
%   0.6 s on (model too fast), because one time constant cannot give two phases.
%   Trains: the first peak is 1.0 where the measurement gives 0.75 to 0.93, and
%   the last peak is too high at 0.5 and 1 Hz; there is no slow state.
%   The sustained current at -80 mV is below the measured 0.096.
%
% OUTPUT:
%   results/fit_report.txt, results/fitted_parameters.mat
%   The frozen values in model/ChR2_Piezo1_opts.m are the restart result,
%   rounded to four significant digits (tau_act 3.162 ms, k_inact 0.9098,
%   tau_rec 0.2023 s). reference_outputs/fit_report.txt is the report of the
%   original fit. Running this script is optional and takes several minutes; it
%   does not change any file of the model.
%===============================================================================

function fit_Piezo1()

here   = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','model'));
outdir = fullfile(here,'..','results');
if ~isfolder(outdir)
    mkdir(outdir);
end
D    = Piezo1_targets();

%% Fit in the logarithm of the three parameters, from a start near the data
p0   = log([4e-3 0.8 0.2]);                    % tau_act (s), k_inact, tau_rec (s)
obj  = @(p) cost(p,D);
opts = optimset('TolX',1e-5,'TolFun',1e-8,'MaxFunEvals',600,'Display','iter');
[p,c]  = fminsearch(obj,p0,opts);
% a second start, to see that the result does not depend on the first
[p2,c2] = fminsearch(obj,log([2e-3 1.0 0.4]),optimset(opts,'Display','off'));

X  = pack(p);
M  = measure_Piezo1(X,D,false);
z  = (M.decay_tau - D.decay.tau)./D.decay.err;

%% Report
fid = fopen(fullfile(outdir,'fit_report.txt'),'w');
say(fid,'Fit of the modified Piezo1 model, written by fit_Piezo1.m on %s\n',datestr(now,'yyyy-mm-dd HH:MM'));
say(fid,'Recording temperature %g K. Cost %.4f (second start %.4f).\n\n',D.T_rec,c,c2);
say(fid,'                 fitted      second start\n');
say(fid,'tau_act (ms)   %8.4f    %8.4f\n',1e3*exp(p(1)),1e3*exp(p2(1)));
say(fid,'k_inact        %8.4f    %8.4f\n',exp(p(2)),exp(p2(2)));
say(fid,'tau_rec (s)    %8.4f    %8.4f\n\n',exp(p(3)),exp(p2(3)));
say(fid,'FITTED TARGETS\n');
say(fid,'%-22s %8s %10s %8s %8s %10s   %s\n','data set','V (mV)','measured','error','model','distance','within the error bar');
for i = 1:numel(D.decay.tau)
    say(fid,'%-22s %8.0f %10.2f %8.2f %8.2f %+10.2f   %s\n',D.decay.names{D.decay.set(i)},D.decay.V_mV(i), ...
        D.decay.tau(i),D.decay.err(i),M.decay_tau(i),z(i),yn(abs(z(i)) <= 1));
end
say(fid,'%-22s %8s %10.2f %8.2f %8.2f %+10.2f   %s\n','deactivation (ms)','-80',D.deact.tau,D.deact.tol,M.deact_tau, ...
    (M.deact_tau-D.deact.tau)/D.deact.tol,yn(abs(M.deact_tau-D.deact.tau) <= D.deact.tol));
say(fid,'%-22s %8s %10.3f %8.3f %8.3f %+10.2f   %s\n','half recovery (s)','-80',D.half.delay,D.half.tol,M.half, ...
    (M.half-D.half.delay)/D.half.tol,yn(abs(M.half-D.half.delay) <= D.half.tol));
npass = sum(abs(z) <= 1) + (abs(M.deact_tau-D.deact.tau) <= D.deact.tol) + (abs(M.half-D.half.delay) <= D.half.tol);
say(fid,'\n%d of %d fitted targets pass. "distance" is in units of the error bar or of the 20 %% tolerance.\n',npass,numel(z)+2);

%% Restart from the first result, with a fresh simplex and tighter limits
% If the cost does not fall and the parameters do not move, the result is a
% stable local minimum. A wider search is not made, so this does not show that
% it is the lowest one.
[pr,cr] = fminsearch(obj,p,optimset('TolX',1e-6,'TolFun',1e-9,'MaxFunEvals',300,'Display','off'));
say(fid,'\nRESTART from the first result\n');
say(fid,'first   : tau_act %.4f ms, k_inact %.4f, tau_rec %.4f s, cost %.5f\n',1e3*exp(p(1)),exp(p(2)),exp(p(3)),c);
say(fid,'restart : tau_act %.4f ms, k_inact %.4f, tau_rec %.4f s, cost %.5f\n',1e3*exp(pr(1)),exp(pr(2)),exp(pr(3)),cr);
if cr < c
    p = pr;
    c = cr;
    say(fid,'The restart is better and replaces the first result.\n');
else
    say(fid,'The restart is not better. The first result is kept.\n');
end
fclose(fid);

fitted = struct('tau_act',exp(p(1)),'k_inact',exp(p(2)),'tau_rec',exp(p(3)),'cost',c,'second_start',exp(p2),'date',datestr(now,'yyyy-mm-dd'));
save(fullfile(outdir,'fitted_parameters.mat'),'fitted');
end

%===============================================================================
% Local functions
%===============================================================================

function X = pack(p)
X = ChR2_Piezo1_opts('tau_act',exp(p(1)),'k_inact',exp(p(2)),'tau_rec',exp(p(3)));
end

function c = cost(p, D)
% Sum of squared distances, each in units of its error bar or tolerance.
X = pack(p);
M = measure_Piezo1(X,D,false);
c = sum(((M.decay_tau - D.decay.tau)./D.decay.err).^2) ...
    + ((M.deact_tau - D.deact.tau)/D.deact.tol)^2 ...
    + ((M.half - D.half.delay)/D.half.tol)^2;
end

function s = yn(b)
if b
    s = 'yes';
else
    s = 'no';
end
end

function say(fid, varargin)
fprintf(varargin{:});
fprintf(fid,varargin{:});
end
