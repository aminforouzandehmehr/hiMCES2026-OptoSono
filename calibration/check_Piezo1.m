%===============================================================================
% Checks the frozen modified Piezo1 model against measurements that were not fitted,
% and draws the channel figure.
%
%   M = check_Piezo1()
%   M : the measurements (measure_Piezo1.m with all checks)
%
% The parameters are read from model/ChR2_Piezo1_opts.m and are not changed.
% reference_outputs/check_report.txt is the report of this script for the
% frozen parameters.
% The rule and the expectation for every check are written in fit_Piezo1.m, before
% the fit. Every check stays in the table, passed or not.
%
% The checks are reported by protocol. The values of one protocol are related
% measurements, so their counts are not added up into one total.
% One protocol is partly calibrated: the fitted delay of half recovery was read
% from the recovery curve that is checked here. The value of that curve next to
% it (0.16 s) is written as tied to the fit and is not counted. The shape of
% the rest of the curve remains a check.
%
% OUTPUT:
%   results/check_report.txt, results/Piezo1_channel_checks.png
%===============================================================================

function M = check_Piezo1()

here   = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','model'));
outdir = fullfile(here,'..','results');
if ~isfolder(outdir)
    mkdir(outdir);
end
D    = Piezo1_targets();
X    = ChR2_Piezo1_opts();
if any(isnan([X.tau_act X.k_inact X.tau_rec]))
    error('the gate parameters are not frozen in ChR2_Piezo1_opts.m');
end
M    = measure_Piezo1(X,D,true);

%% Report
fid = fopen(fullfile(outdir,'check_report.txt'),'w');
say(fid,'Checks of the modified Piezo1 model, written by check_Piezo1.m\n');
say(fid,'tau_act %.4f ms, k_inact %.4f, tau_rec %.4f s, at %g K. None of the values below entered the fit.\n',1e3*X.tau_act,X.k_inact,X.tau_rec,D.T_rec);
say(fid,'A check passes within two standard errors or two standard deviations, or within 20 %% where the source gives no error.\n');
say(fid,['One protocol is partly calibrated: the fitted delay of half recovery (%.3f s) was read from the recovery curve\n' ...
         'checked below. The value of that curve next to it (0.16 s) is tied to the fit and is not counted.\n\n'],D.half.delay);
say(fid,'%-52s %10s %10s %10s   %s\n','','measured','allowed','model','passes');

say(fid,'Inactivation of human PIEZO1, whole cell\n');
pHuman = line(fid,'  at -80 mV (ms)',D.human.tau,2*D.human.sem,M.human_tau);

say(fid,'\nInactivation of mouse Piezo1 below the fitted voltages, whole cell\n');
neg    = D.wufar.V_mV < 0;
pNeg   = line(fid,'  at -80 mV, Zheng (ms)',D.zheng80.tau,2*D.zheng80.sem,M.zheng80);
for i = find(neg)
    pNeg(end+1) = line(fid,sprintf('  at %+d mV, Wu (ms)',D.wufar.V_mV(i)),D.wufar.tau(i),D.wufar.tol(i),M.wufar_tau(i)); %#ok<AGROW>
end

say(fid,'\nInactivation of mouse Piezo1 above the fitted voltages, whole cell\n');
pPos   = false(1,0);
for i = find(~neg)
    pPos(end+1) = line(fid,sprintf('  at %+d mV, Wu (ms)',D.wufar.V_mV(i)),D.wufar.tau(i),D.wufar.tol(i),M.wufar_tau(i)); %#ok<AGROW>
end

say(fid,'\nSustained current, cell attached, -80 mV\n');
pSust  = line(fid,'  current left after 300 ms (fraction)',D.sust.value,2*D.sust.sd,M.sust);

say(fid,'\nRecovery after two steps of 100 ms, cell attached, -80 mV (partly calibrated, see above)\n');
tied   = abs(D.rec.delay - 0.16) < 1e-9;
pRec   = false(1,0);
for i = 1:numel(D.rec.delay)
    p = line(fid,sprintf('  recovery after %.2f s',D.rec.delay(i)),D.rec.value(i),2*D.rec.sem(i),M.rec(i),tied(i));
    if ~tied(i)
        pRec(end+1) = p; %#ok<AGROW>
    end
end

say(fid,'\nTrains of 30 ms steps for 4 s, as a fraction of a 300 ms step given 10 s earlier, cell attached\n');
pFirst = false(1,0);
pLast  = false(1,0);
for i = 1:numel(D.train.f_Hz)
    pFirst(end+1) = line(fid,sprintf('  first peak at %g Hz',D.train.f_Hz(i)),D.train.first(i),2*D.train.first_sem(i),M.train_first(i)); %#ok<AGROW>
end
for i = 1:numel(D.train.f_Hz)
    pLast(end+1) = line(fid,sprintf('  last peak at %g Hz',D.train.f_Hz(i)),D.train.last(i),2*D.train.last_sem(i),M.train_last(i)); %#ok<AGROW>
end

say(fid,'\nBY PROTOCOL\n');
say(fid,'The values of one protocol are related measurements. The counts are given by protocol and are not added up.\n');
say(fid,'%-84s %s\n','human PIEZO1, inactivation at -80 mV (Albuisson et al. 2013)',cnt(pHuman));
say(fid,'%-84s %s\n','mouse Piezo1, inactivation at -100 and -80 mV (Wu et al. 2017, Zheng et al. 2019)',cnt(pNeg));
say(fid,'%-84s %s\n','mouse Piezo1, inactivation at +40 to +100 mV (Wu et al. 2017)',cnt(pPos));
say(fid,'%-84s %s\n','sustained current after 300 ms (Lewis et al. 2017)',cnt(pSust));
say(fid,'%-84s %s\n','recovery curve, values not tied to the fit (Lewis et al. 2017)',cnt(pRec));
say(fid,'%-84s %s\n','trains, first peak (Lewis et al. 2017)',cnt(pFirst));
say(fid,'%-84s %s\n','trains, last peak (Lewis et al. 2017)',cnt(pLast));
say(fid,'\nOpen fraction at the peak of a saturating poke at -80 mV: %.3f.\n',M.peak_open);
say(fid,['Not represented in this form, and not run: a slow component of recovery, and the\n' ...
         'recovery of the current by an outward current at positive potentials.\n']);
fclose(fid);

%% Figure
ink  = [0.04 0.04 0.04];
sub  = [0.32 0.32 0.31];
blue = [42 120 214]/255;
ora  = [235 104 52]/255;
grey = [0.62 0.62 0.61];
f    = figure('Color','w','Units','inches','Position',[0 0 11 8.5],'Visible','off');
tl   = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
sty  = @(ax) set(ax,'FontSize',9.5,'TickDir','out','Box','off','LineWidth',0.6);

% a: decay time constant against voltage
ax = nexttile(tl);
hold(ax,'on');
sty(ax);
Vf   = -100:10:100;
tauB = arrayfun(@(v) one_decay(X,v,D.T_rec),Vf);
h3   = plot(ax,Vf,tauB,'-','Color',blue,'LineWidth',1.6);
k1   = D.decay.set == 1;
h1   = errorbar(ax,D.decay.V_mV(k1)-1.5,D.decay.tau(k1),D.decay.err(k1),'o','Color',ink,'MarkerFaceColor',ink,'MarkerSize',5,'CapSize',3);
h2   = errorbar(ax,D.decay.V_mV(~k1)+1.5,D.decay.tau(~k1),D.decay.err(~k1),'s','Color',ora,'MarkerFaceColor',ora,'MarkerSize',5,'CapSize',3);
h4   = plot(ax,D.wufar.V_mV,D.wufar.tau,'o','Color',ink,'MarkerSize',5);
h5   = errorbar(ax,-80,D.human.tau,2*D.human.sem,'d','Color',grey,'MarkerFaceColor',grey,'MarkerSize',6,'CapSize',3);
set(ax,'YScale','log');
ylim(ax,[5 120]);
yticks(ax,[5 10 20 50 100]);
xlabel(ax,'V (mV)');
ylabel(ax,'decay time constant of the current (ms)');
title(ax,'Inactivation against voltage','FontSize',11,'FontWeight','bold');
legend(ax,[h1 h2 h4 h5 h3],{'Wu et al. 2017, fitted','Zheng et al. 2019, fitted','Wu et al. 2017, checked only','human, Albuisson et al. 2013, checked only','model'}, ...
    'Location','northwest','Box','off','FontSize',8.5);

% b: a held stimulus at -80 mV, the modified model against version 1.1
ax = nexttile(tl);
hold(ax,'on');
sty(ax);
t   = (0:2e-5:0.3)';
S   = channel_run(X,-80,D.T_rec,[0 0.001 0.3],[0 2*X.F_ref 2*X.F_ref],t);
rel = 1 - exp(-t/0.1);                        % released code: the activation gate alone, tau 0.1 s
h1  = plot(ax,1e3*t,rel/max(rel),'-','Color',grey,'LineWidth',1.6);
h2  = plot(ax,1e3*t,S.g/max(S.g),'-','Color',blue,'LineWidth',1.6);
xlabel(ax,'time after stimulus onset (ms)');
ylabel(ax,'current (fraction of its largest value)');
ylim(ax,[0 1.05]);
title(ax,'A held stimulus at -80 mV','FontSize',11,'FontWeight','bold');
legend(ax,[h1 h2],{'version 1.1','modified model'},'Location','east','Box','off','FontSize',8.5);

% c: recovery
ax = nexttile(tl);
hold(ax,'on');
sty(ax);
df  = logspace(-2,log10(41),60);
rf  = arrayfun(@(d) one_recovery(X,D.T_rec,d),df);
h2  = plot(ax,df,rf,'-','Color',blue,'LineWidth',1.6);
h1  = errorbar(ax,D.rec.delay,D.rec.value,D.rec.sem,'o','Color',ink,'MarkerFaceColor',ink,'MarkerSize',5,'CapSize',3);
h3  = plot(ax,D.half.delay,0.5,'x','Color',ora,'MarkerSize',10,'LineWidth',1.8);
set(ax,'XScale','log');
xlim(ax,[0.008 50]);
ylim(ax,[0 1.05]);
xlabel(ax,'delay between the two steps (s)');
ylabel(ax,'recovery');
title(ax,'Recovery from inactivation','FontSize',11,'FontWeight','bold');
legend(ax,[h1 h3 h2],{'Lewis et al. 2017, shape of the curve checked','half recovery, fitted, read from this curve','model'},'Location','southeast','Box','off','FontSize',8.5);

% d: trains
ax = nexttile(tl);
hold(ax,'on');
sty(ax);
h1 = errorbar(ax,D.train.f_Hz,D.train.first,D.train.first_sem,'o','Color',ink,'MarkerSize',5,'CapSize',3);
h2 = errorbar(ax,D.train.f_Hz,D.train.last,D.train.last_sem,'o','Color',ink,'MarkerFaceColor',ink,'MarkerSize',5,'CapSize',3);
h3 = plot(ax,D.train.f_Hz,M.train_first,'--','Color',blue,'LineWidth',1.4);
h4 = plot(ax,D.train.f_Hz,M.train_last,'-','Color',blue,'LineWidth',1.6);
set(ax,'XScale','log');
xlim(ax,[0.4 12]);
xticks(ax,D.train.f_Hz);
ylim(ax,[0 1.1]);
xlabel(ax,'rate of the 30 ms steps (Hz)');
ylabel(ax,'peak (fraction of a single step)');
title(ax,'Repeated stimuli for 4 s','FontSize',11,'FontWeight','bold');
legend(ax,[h1 h2 h3 h4],{'measured, first peak','measured, last peak','model, first peak','model, last peak'},'Location','southwest','Box','off','FontSize',8.5);

title(tl,'Modified Piezo1 model, the isolated channel against measurements','FontWeight','bold','FontSize',13);
subtitle(tl,'room temperature; inactivation from whole cell recordings, recovery and trains from cell attached recordings (Lewis et al. 2017); error bars as plotted in the sources', ...
    'FontSize',9.5,'Color',sub);
exportgraphics(f,fullfile(outdir,'Piezo1_channel_checks.png'),'Resolution',220);
close(f);
end

%===============================================================================
% Local functions
%===============================================================================

function pass = line(fid, label, meas, allowed, model, tied)
% Writes one check and returns whether it passes. A value tied to the fit is
% written as such.
pass = abs(model-meas) <= allowed;
if nargin > 5 && tied
    s = 'tied to the fit, not counted';
elseif pass
    s = 'yes';
else
    s = 'no';
end
say(fid,'%-52s %10.3f %10.3f %10.3f   %s\n',label,meas,allowed,model,s);
end

function s = cnt(p)
% Number of passed values of one protocol.
s = sprintf('%d of %d',sum(p),numel(p));
end

function tau = one_decay(X, V_mV, T)
% Decay time constant in ms from the peak of a 300 ms poke, as in measure_Piezo1.m.
t      = (0:2e-5:0.301)';
S      = channel_run(X,V_mV,T,[0 0.001 0.301],[0 2*X.F_ref 2*X.F_ref],t);
[~,ip] = max(S.g);
tt     = S.t(ip:end)-S.t(ip);
yy     = S.g(ip:end);
tau    = 1e3*fminbnd(@(x) res(x,tt,yy),2e-4,2,optimset('TolX',1e-9));
end

function s = res(tau, t, y)
A = [ones(size(t)) exp(-t/tau)];
c = A\y;
s = sum((y - A*c).^2);
end

function r = one_recovery(X, T, delay)
% Recovery after two steps of 100 ms, as in measure_Piezo1.m.
ramp = 0.002;
d    = 0.100;
F    = 2*X.F_ref;
t2   = d + ramp + delay;
tb   = [0 ramp d d+ramp t2 t2+ramp t2+d];
Fb   = [0 F F 0 0 F F];
tout = [(0:2e-5:d)'; (t2:2e-5:t2+d)'];
S    = channel_run(X,-80,T,tb,Fb,tout);
w1   = S.t <= d;
p1   = max(S.g(w1));
e1   = S.g(find(w1,1,'last'));
p2   = max(S.g(~w1));
r    = (p2 - e1)/(p1 - e1);
end

function say(fid, varargin)
fprintf(varargin{:});
fprintf(fid,varargin{:});
end
