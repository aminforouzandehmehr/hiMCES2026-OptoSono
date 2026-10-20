%===============================================================================
% Settings of the two simulated cells of Fig. 8A
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
%   X = OptoSono_settings(cellType)
%   cellType : 'Piezo1'     hiMCES2026 with ChR2 and the modified Piezo1 model (40 states)
%              'Piezo1off'  the same cell with the Piezo1 current off (39 states)
%
% Both cells: Tyrode solution (classicVSoptic = 1) at 304 K (31 C); light pulses
% of 10 ms at 0.15 mW/mm2 and 0.5 Hz; applied mechanical pulses of 0.5 s. The
% number and timing of the light and mechanical pulses are set by the caller.
%
% Piezo1 cell: 300 times the published conductance (91.17 S/F, about 9 nS per
% cell), reference stress F_ref 0.336 kPa (half activation at 0.168 kPa), native
% stress filter of 0.15 s, gate kinetics as fitted (ChR2_Piezo1_opts.m), 5 % of
% the inward Piezo1 current as Ca2+.
%===============================================================================

function X = OptoSono_settings(cellType)

switch cellType
    case 'Piezo1'
        X = ChR2_Piezo1_opts('T_override',304,'gp_scale',300,'F_ref',0.336,'tau_trans',0.15);
    case 'Piezo1off'
        X = ChR2_Piezo1_opts('T_override',304,'gp_scale',0);
    otherwise
        error('cellType must be ''Piezo1'' or ''Piezo1off''');
end

%% Light pulses (ChR2)
X.optIrr    = 0.15;       % mW/mm2
X.optDur    = 0.010;      % second
X.optFreq   = 0.5;        % hertz
X.optT0     = 0;          % second, set by the caller
X.optN      = 0;          % set by the caller

%% Applied mechanical pulses (Piezo1)
X.touchF    = 0;          % kPa, set by the caller
X.touchT0   = 0;          % second, set by the caller
X.touchDur  = 0.5;        % second
X.touchFreq = 0.5;        % hertz
X.touchN    = 0;          % set by the caller
end
