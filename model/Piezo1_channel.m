%===============================================================================
% The two gates of the modified Piezo1 model on their own
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
%   [dma, dh, r] = Piezo1_channel(V, ma, h, F_Pie, T, X)
%
%   V      : membrane potential in volt
%   ma     : activation gate (0 to 1)
%   h      : availability (0 to 1); 1 means no channel is inactivated
%   F_Pie  : stress that the activation gate senses, kPa
%   T      : temperature in kelvin
%   X      : options from ChR2_Piezo1_opts
%   dma, dh: the two derivatives per second
%   r      : steady state of the activation gate and the three time constants
%
% These are the statements of the Piezo1 block of hiMCES2026_ChR2_Piezo1.m for
% the modified gating (X.gating = 'reduced'). verify_release.m checks that the
% two agree. The calibration scripts use this function, so they need no cell
% model.
%
% Modified from Piezo1.m version 1.1
% (https://github.com/aminforouzandehmehr/Piezo1-Channel-Model) of
% Forouzandehmehr et al., 2023 (https://doi.org/10.22489/CinC.2023.359).
%===============================================================================

function [dma, dh, r] = Piezo1_channel(V, ma, h, F_Pie, T, X)

q_Pie  = 2.0^((X.T_rec - T)/10.0);   % dimensionless, 1 at the recording temperature
d      = (F_Pie/X.F_ref) * 5;        % dimensionless
mpi    = 1.0/(1+exp(-(d-2.5)/0.25));
al_Xs  = 1/sqrt(1.0+exp((-60.0+V*1000.0)/50.0));
be_Xs  = 900.0/(1.0+exp((+171.0-V*1000.0)/72.0));
tau_p  = X.tau_act*q_Pie;                                % second, activation and deactivation
ta_Xs  = 1e-3 * X.k_inact*q_Pie * 30.0*al_Xs*be_Xs/100.0;   % second, inactivation of an open channel
tau_r  = X.tau_rec*q_Pie;                                % second, recovery from inactivation

dma    = (mpi-ma)/tau_p;
dh     = (1.0-h)/tau_r - ma*h/ta_Xs;

r.ma_inf    = mpi;
r.tau_act   = tau_p;
r.tau_inact = ta_Xs;
r.tau_rec   = tau_r;
end
