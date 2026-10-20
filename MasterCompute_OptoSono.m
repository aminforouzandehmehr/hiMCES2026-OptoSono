%===============================================================================
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
% PURPOSE:
%   Reproduces Fig. 8A of Pettebone, Han et al., bioRxiv 2026
%   (DOI: 10.64898/2026.10.19.760310) with hiMCES2026, ChR2 and the modified
%   Piezo1 model: three light pulses, then one applied mechanical pulse.
%   Controls: no mechanical pulse; Piezo1 current off.
%
% MAIN FLAGS:
%   stimAmp   : effective stress at the Piezo1 gate, kPa (Fig. 8A: 0.336)
%   stimOn    : onset of the mechanical pulse, s (Fig. 8A: 5.0)
%   stimDur   : duration of the mechanical pulse, s (Fig. 8A: 0.5)
%   saveFigs  : 1 writes results/Figure8A.pdf and results/Figure8A.png
%
% INITIAL CONDITIONS:
%   initial_conditions/OptoSono_31C_Piezo1.mat     40 states
%   initial_conditions/OptoSono_31C_Piezo1off.mat  39 states
%   Both are states at the onset of a light pulse after 1000 s of 0.5 Hz
%   light pacing at 304 K (condition_OptoSono.m).
%===============================================================================

%% Housekeeping
clear; close all; clc
tic
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'model'));
addpath(fullfile(here,'analysis'));

%% ============================================================================
%  1) FLAGS
% ============================================================================
stimAmp  = 0.336;       % kPa
stimOn   = 5.0;         % second, 1 s after the third light pulse
stimDur  = 0.5;         % second
saveFigs = 0;           % 0 no, 1 yes

%% ============================================================================
%  2) INITIAL CONDITIONS
% ============================================================================
IC1 = load(fullfile(here,'initial_conditions','OptoSono_31C_Piezo1.mat'));
IC0 = load(fullfile(here,'initial_conditions','OptoSono_31C_Piezo1off.mat'));

%% ============================================================================
%  3) PARAMETERS
% ============================================================================
classicVSoptic = 1;     % Tyrode solution; temperature 304 K is set in OptoSono_settings.m
dt             = 2e-4;  % second, output interval
tLight         = 4.9;   % second, end of the light part (pulses at 0, 2 and 4 s)
tEnd           = 9.5;   % second

if stimOn < tLight
    error('stimOn must be at or after %.1f s',tLight);
end

X1 = OptoSono_settings('Piezo1');
X0 = OptoSono_settings('Piezo1off');

% Light part: three pulses from t = 0
X1L      = X1;
X1L.optN = 3;
X0L      = X0;
X0L.optN = 3;

% Mechanical part: one pulse, no light
X1M          = X1;
X1M.touchF   = stimAmp;
X1M.touchT0  = stimOn;
X1M.touchDur = stimDur;
X1M.touchN   = 1;
X0M          = X0;
X0M.touchF   = stimAmp;
X0M.touchT0  = stimOn;
X0M.touchDur = stimDur;
X0M.touchN   = 1;
X1N          = X1;      % control without the mechanical pulse

%% ============================================================================
%  4) INTEGRATION
% ============================================================================
light1 = run_segments(X1L, IC1.Y,       0,      tLight, dt, classicVSoptic);
stim1  = run_segments(X1M, light1.Yend, tLight, tEnd,   dt, classicVSoptic);
nost1  = run_segments(X1N, light1.Yend, tLight, tEnd,   dt, classicVSoptic);
light0 = run_segments(X0L, IC0.Y,       0,      tLight, dt, classicVSoptic);
stim0  = run_segments(X0M, light0.Yend, tLight, tEnd,   dt, classicVSoptic);
toc

%% ============================================================================
%  5) POST-PROCESS
% ============================================================================
M = figure8A_metrics(light1, stim1, nost1, stim0, stimOn, stimDur);

%% ============================================================================
%  6) PLOTS
% ============================================================================
BLUE = [0.10 0.45 0.85];
ORNG = [0.90 0.42 0.10];
GRN  = [0.15 0.55 0.30];
GRY  = [0.60 0.60 0.60];
fs   = 8;
lw   = 1.2;

t  = [light1.t; stim1.t(stim1.t > light1.t(end))];
V  = [light1.Vm_mV; stim1.Vm_mV(stim1.t > light1.t(end))];
C  = [light1.Cai_nM; stim1.Cai_nM(stim1.t > light1.t(end))];

fig = figure('Color','w','Units','inches','Position',[1 1 7.2 6.0]);
tl  = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');

% a: membrane potential
ax = nexttile(tl,[1 3]);
hold(ax,'on');
yl = [-90 65];
for p = [0 2 4]
    patch(ax,[p p p+X1.optDur p+X1.optDur],[yl(1) yl(2) yl(2) yl(1)],BLUE,'EdgeColor',BLUE,'LineWidth',1.5);
end
patch(ax,[stimOn stimOn stimOn+stimDur stimOn+stimDur],[yl(1) yl(2) yl(2) yl(1)],ORNG,'FaceAlpha',0.25,'EdgeColor','none');
plot(ax,nost1.t,nost1.Vm_mV,'--','Color',GRY,'LineWidth',1.0);
plot(ax,t,V,'k','LineWidth',lw);
xlim(ax,[-0.3 tEnd]);
ylim(ax,yl);
xlabel(ax,'time (s)');
ylabel(ax,'V_m (mV)');
set(ax,'box','off','tickdir','out','fontsize',fs);

% b: calcium
ax = nexttile(tl,[1 3]);
hold(ax,'on');
plot(ax,nost1.t,nost1.Cai_nM,'--','Color',GRY,'LineWidth',1.0);
plot(ax,t,C,'Color',GRN,'LineWidth',lw);
xlim(ax,[-0.3 tEnd]);
xlabel(ax,'time (s)');
ylabel(ax,'[Ca^{2+}]_i (nM)');
set(ax,'box','off','tickdir','out','fontsize',fs);

% c: ChR2 current of the third light pulse
ax = nexttile(tl);
w  = light1.t >= 3.99 & light1.t <= 4.20;
plot(ax,1e3*(light1.t(w)-4),light1.I_ChR2(w),'Color',BLUE,'LineWidth',lw);
xlim(ax,[-10 200]);
xlabel(ax,'time from light onset (ms)');
ylabel(ax,'I_{ChR2} (A/F)');
set(ax,'box','off','tickdir','out','fontsize',fs);

% d: Piezo1 current of the mechanical pulse, with the first 30 ms enlarged
ax = nexttile(tl);
w  = stim1.t >= stimOn-0.01 & stim1.t <= stimOn+stimDur+0.1;
plot(ax,1e3*(stim1.t(w)-stimOn),stim1.I_Piezo1(w),'Color',ORNG,'LineWidth',lw);
xlim(ax,[-10 1e3*stimDur+100]);
xlabel(ax,'time from stimulus onset (ms)');
ylabel(ax,'I_{Piezo1} (A/F)');
set(ax,'box','off','tickdir','out','fontsize',fs);
drawnow
pos = ax.Position;
ai  = axes(fig,'Position',[pos(1)+0.50*pos(3) pos(2)+0.25*pos(4) 0.45*pos(3) 0.40*pos(4)]);
w   = stim1.t >= stimOn & stim1.t <= stimOn+0.03;
plot(ai,1e3*(stim1.t(w)-stimOn),stim1.I_Piezo1(w),'Color',ORNG,'LineWidth',lw);
xlim(ai,[0 30]);
set(ai,'box','on','fontsize',fs-1);

% e: light and mechanical beats aligned at the upstroke
ax = nexttile(tl);
hold(ax,'on');
ul = upstroke_times(light1.t,light1.Vm_mV);
us = upstroke_times(stim1.t,stim1.Vm_mV);
us = us(us >= stimOn);
w  = light1.t >= ul(2)-0.05 & light1.t <= min(ul(2)+1.85,light1.t(end));
plot(ax,1e3*(light1.t(w)-ul(2)),light1.Vm_mV(w),'Color',BLUE,'LineWidth',lw);
if ~isempty(us)
    w = stim1.t >= us(1)-0.05 & stim1.t <= us(1)+1.85;
    plot(ax,1e3*(stim1.t(w)-us(1)),stim1.Vm_mV(w),'Color',ORNG,'LineWidth',lw);
end
xlim(ax,[-50 1850]);
ylim(ax,yl);
xlabel(ax,'time from upstroke (ms)');
ylabel(ax,'V_m (mV)');
set(ax,'box','off','tickdir','out','fontsize',fs);

if saveFigs == 1
    outdir = fullfile(here,'results');
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    exportgraphics(fig,fullfile(outdir,'Figure8A.pdf'),'ContentType','vector');
    exportgraphics(fig,fullfile(outdir,'Figure8A.png'),'Resolution',300);
end
