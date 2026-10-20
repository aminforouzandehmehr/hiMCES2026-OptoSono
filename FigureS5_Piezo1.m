%===============================================================================
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
% PURPOSE:
%   Reproduces Supplementary Fig. S5 of Pettebone, Han et al., bioRxiv 2026
%   (DOI: 10.64898/2026.10.19.760310): decay time constant of the Piezo1 current
%   against voltage for the isolated channel of the modified Piezo1 model at the
%   recording temperature (295 K), with the data it was fitted to and the data
%   used only as checks. Prints the fitted targets and the checks.
%   No parameter is fitted here.
%
% MAIN FLAGS:
%   saveFigs  : 1 writes results/FigureS5.pdf and results/FigureS5.png
%===============================================================================

%% Housekeeping
clear; close all; clc
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'model'));
addpath(fullfile(here,'calibration'));

%% ============================================================================
%  1) FLAGS
% ============================================================================
saveFigs = 0;           % 0 no, 1 yes

%% ============================================================================
%  2) PARAMETERS AND DATA
% ============================================================================
X  = ChR2_Piezo1_opts();            % frozen gate parameters of the modified model
D  = Piezo1_targets();              % measured values, with sources
Dc = D;
Dc.decay.V_mV = -100:5:100;         % mV, dense curve

%% ============================================================================
%  3) CHANNEL PROTOCOLS
% ============================================================================
Mc = measure_Piezo1(X,Dc,false);    % curve
M  = measure_Piezo1(X,D,true);      % fitted targets and checks

%% ============================================================================
%  4) REPORT
% ============================================================================
fprintf('Gate parameters at %d K: tau_act %.3f ms, k_inact %.4f, tau_rec %.4f s\n', ...
        D.T_rec,1e3*X.tau_act,X.k_inact,X.tau_rec);
fprintf('Fitted:  decay at -60, -40, -20, 0, +20 mV, model %s ms\n',mat2str(round(M.decay_tau(1:5),2)));
fprintf('         deactivation at -80 mV, measured %.1f ms, model %.2f ms\n',D.deact.tau,M.deact_tau);
fprintf('         half recovery at -80 mV, measured %.3f s, model %.3f s\n',D.half.delay,M.half);
fprintf('Checks:  human PIEZO1 at -80 mV, measured %.1f +/- %.1f ms, model %.2f ms\n',D.human.tau,D.human.sem,M.human_tau);
fprintf('         mouse at -80 mV (Zheng), measured %.1f ms, model %.2f ms\n',D.zheng80.tau,M.zheng80);
fprintf('         sustained after 300 ms, measured %.3f (s.d. %.3f), model %.3f\n',D.sust.value,D.sust.sd,M.sust);
fprintf('         decay +40 to +100 mV, measured %s ms, model %s ms\n',mat2str(D.wufar.tau(3:6)),mat2str(round(M.wufar_tau(3:6),1)));
fprintf('The full table of checks, including recovery and repeated stimuli, is written by calibration/check_Piezo1.m.\n');

%% ============================================================================
%  5) PLOT
% ============================================================================
BLK = [0 0 0];
ORG = [0.85 0.40 0.10];
GRY = [0.55 0.55 0.55];
TEA = [8 126 128]/255;
fs  = 8;

fig = figure('Color','w','Units','inches','Position',[1 1 3.6 3.0]);
ax  = axes(fig);
hold(ax,'on');
w1 = D.decay.set == 1;
w2 = D.decay.set == 2;
h1 = errorbar(ax,D.decay.V_mV(w1)-1.5,D.decay.tau(w1),D.decay.err(w1),'o','Color',BLK, ...
              'MarkerFaceColor',BLK,'MarkerSize',4,'CapSize',2);
h2 = errorbar(ax,D.decay.V_mV(w2)+1.5,D.decay.tau(w2),D.decay.err(w2),'s','Color',ORG, ...
              'MarkerFaceColor',ORG,'MarkerSize',4,'CapSize',2);
h3 = plot(ax,D.wufar.V_mV,D.wufar.tau,'o','Color',BLK,'MarkerFaceColor','w','MarkerSize',4);
h4 = errorbar(ax,-80,D.human.tau,2*D.human.sem,'d','Color',GRY,'MarkerFaceColor',GRY,'MarkerSize',5,'CapSize',2);
h5 = plot(ax,Mc.decay_V,Mc.decay_tau,'Color',TEA,'LineWidth',1.4);
set(ax,'YScale','log','YTick',[5 10 20 50 100],'box','off','tickdir','out','fontsize',fs);
xlim(ax,[-105 105]);
ylim(ax,[5 120]);
xlabel(ax,'V (mV)');
ylabel(ax,'decay time constant (ms)');
legend(ax,[h1 h2 h3 h4 h5],{'Wu et al. 2017, fit','Zheng et al. 2019, fit','Wu et al. 2017, check', ...
       'Albuisson et al. 2013 (human), check','model'},'Location','northwest','Box','off','FontSize',fs-1);

if saveFigs == 1
    outdir = fullfile(here,'results');
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    exportgraphics(fig,fullfile(outdir,'FigureS5.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outdir,'FigureS5.png'),'Resolution',300);
end
