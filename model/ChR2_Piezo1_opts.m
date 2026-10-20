%===============================================================================
% Options and parameters of hiMCES2026_ChR2_Piezo1
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
% PURPOSE:
%   Returns the structure X that the model takes as its last argument. Any field
%   can be changed with a name and value pair, for example
%       X = ChR2_Piezo1_opts('T_override', 304, 'gp_scale', 100);
%
%   With gating = 'released' and accounting = 'published' the Piezo1 block is
%   Piezo1.m version 1.1. With gp_scale = 0 and no light the model is the
%   published hiMCES2026.
%
%   The defaults below are not the configuration of Fig. 8A: conductance and
%   F_ref are the published values of version 1.1, and no stress filter is set.
%   OptoSono_settings.m returns the configuration of Fig. 8A.
%
% MAIN FIELDS:
%   T_override      : temperature in kelvin, [] keeps the published choice
%   gp_scale        : Piezo1 conductance relative to the published 0.3039 S/F
%   gating          : 'reduced' modified model, 'released' version 1.1
%   accounting      : 'split' between Na+ and Ca2+ (simplified, see below),
%                     'published' all as Ca2+
%   F_ref           : reference stress of the activation gate, kPa (half
%                     activation at F_ref/2)
%   tau_trans       : time constant of the native stress filter, s ([] = none)
%   touchF          : applied mechanical stimulus, an effective stress added to
%                     the stress that the Piezo1 activation gate senses, kPa
%                     (0 = none). It is not an acoustic pressure.
%   optIrr          : irradiance during a light pulse, mW/mm2 (0 = dark)
%
% STATUS:
%   The three gate parameters of the modified model were fitted to published
%   channel recordings by calibration/fit_Piezo1.m and are frozen below.
%===============================================================================

function X = ChR2_Piezo1_opts(varargin)

%% Temperature
% The bath is chosen by classicVSoptic in the model call, not here.
X.T_override    = [];        % kelvin; [] keeps the published 310 K or 294 K; 304 K is 31 C

%% Piezo1
X.gp_scale      = 1;         % 1 is the published conductance, 0 switches the current off
X.gating        = 'reduced'; % 'reduced' or 'released'
X.accounting    = 'split';   % 'split' or 'published'
X.F_ref         = 0.056;     % kPa, reference stress of the activation gate, as published

%% Native stress filter, an assumption
% With tau_trans the cell's own stress reaches the activation gate filtered:
%   d(frc_smooth)/dt = (frc - frc_smooth)/tau_trans
% and the gate senses frc_smooth plus the applied stimulus. The applied stimulus
% is not filtered. The model then has a 40th state, frc_smooth. The time constant
% is not measured. An empty value means no filter and 39 states.
X.tau_trans     = [];        % second, [] = the gate senses the stress as it is

%% Modified gating, frozen on 11 October 2026 (calibration/fit_Piezo1.m)
% The second gate is the availability h. It is 1 when no channel is inactivated.
% An open channel inactivates with the time constant tau_inact(V), and
% inactivated channels recover with tau_rec:
%   dh/dt = (1 - h)/tau_rec - ma*h/tau_inact(V)
% tau_inact(V) is k_inact times the published voltage function of the gate time
% constant, in milliseconds. The three values hold at the recording temperature
% T_rec. At another temperature all three are scaled with the Q10 of 2 that the
% published model uses for its own time constants; for Piezo1 this Q10 is an
% assumption.
X.T_rec         = 295;       % kelvin, room temperature of the recordings, taken as 22 C
X.tau_act       = 3.162e-3;  % second, fitted to the measured deactivation; activation assumed equal
X.k_inact       = 0.9098;    % dimensionless, scale on the published voltage function
X.tau_rec       = 0.2023;    % second, recovery from inactivation

%% Accounting of the Piezo1 current, used with accounting = 'split'
% Simplified ion accounting. The net current is divided by its sign, not into
% Na+, K+ and Ca2+ currents with their own driving forces.
% Inward net current is divided between Na+ and Ca2+. Outward net current is K+
% leaving, which enters no balance because the published model keeps K_i constant.
% Ca2+ that would still enter while the net current is outward is not counted.
% f_Ca is the Ca2+ share of the inward current. 0.05 is the share of constant
% field currents at -75 mV in the two baths of the model (0.046 to 0.055), with
% a permeability of Ca2+ relative to Na+ of 1.3. That ratio combines two
% measurements (Ca2+ to Cs+ 1.21 in mouse Piezo1, Coste et al. 2015, Nat Commun
% 6:7223; Cs+ to Na+ 1.07 in human PIEZO1, Gnanasambandam et al. 2015, PLoS One
% 10:e0125503) and is an assumption for this model.
X.f_Ca          = 0.05;      % dimensionless

%% Applied mechanical stimulus
% While a pulse is on, touchF is added to the stress of the myofilament model in
% the input of the activation gate.
X.touchF        = 0;         % kPa, 0 = no stimulus
X.touchT0       = 0;         % second, start of the first pulse
X.touchDur      = 0.050;     % second, duration of each pulse
X.touchFreq     = 0.5;       % hertz
X.touchN        = 0;         % number of pulses

%% Light train (ChR2, Williams et al. 2013)
X.optIrr        = 0;         % mW/mm2, irradiance during a pulse (0 = dark)
X.optT0         = 0;         % second, start of the first pulse
X.optDur        = 0.050;     % second, duration of each pulse
X.optFreq       = 0.5;       % hertz
X.optN          = 0;         % number of pulses
X.chr2_gmax     = 1.0;       % scale on the published ChR2 conductance

%% Name and value pairs given by the caller
for i = 1:2:numel(varargin)
    if ~isfield(X,varargin{i})
        error('unknown option ''%s''',varargin{i});
    end
    X.(varargin{i}) = varargin{i+1};
end
end
