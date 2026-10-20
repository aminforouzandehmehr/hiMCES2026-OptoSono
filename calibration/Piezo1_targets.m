%===============================================================================
% The measured channel responses used for the modified Piezo1 model: three that are
% fitted and six that are only checked after the parameters are frozen.
%
%   D = Piezo1_targets()
%
% All recordings are of Piezo1 expressed in HEK293 cells at room temperature,
% taken here as 22 C (295 K). Two preparations are mixed, because neither covers
% everything: whole cell with a blunt probe (the stimulus is too slow to resolve
% activation) and cell attached with a fast pressure clamp (inactivation is
% about three times slower there than in whole cell).
%
% FITTED:
%   D.decay   inactivation time constant against voltage, mouse Piezo1, whole
%             cell, blunt probe. Wu et al. 2017, Cell Rep 21:2357, Fig 1B, and
%             Zheng et al. 2019, eLife 8:e44003, Fig 2H. The values are the
%             points plotted in Fig 2D of the CinC 2023 paper, read from the
%             paper rendered at 1100 dpi (reading error about 0.1 ms). The Wu
%             points were compared with Fig 1B of that paper and agree within
%             the reading error of that figure. The error bars are as plotted.
%   D.deact   time constant of the current decay after the pressure is released
%             to +5 mmHg, mouse Piezo1, cell attached, -80 mV. Lewis et al.
%             2017, Cell Rep 19:2572, Fig S4I: the fast component, 4.0 ms, which
%             carries 82 % of the decay (Fig S4H). Read from the figure; no
%             error bar is drawn.
%   D.half    delay at which the two step recovery reaches one half. Same paper,
%             Fig S4F: 0.439 at 0.08 s and 0.514 at 0.16 s, so 0.140 s by
%             interpolation on the logarithmic axis. The measured curve is flat
%             there (0.42 to 0.55 between 0.04 and 0.64 s), so this delay is not
%             sharply defined.
%
% CHECKED ONLY (none of them enters the fit; all were inspected during development):
%   D.human   inactivation time constant of human PIEZO1, whole cell, -80 mV:
%             8.6 +/- 0.4 ms (n = 25), Albuisson et al. 2013, Nat Commun 4:1884.
%   D.zheng80 mouse Piezo1 at -80 mV: 11.9 +/- 0.6 ms, Zheng et al. 2019, text.
%   D.wufar   Wu et al. 2017, Fig 1B outside the fitted range, read by eye from
%             the figure at low resolution (about 10 %).
%   D.rec     the whole two step recovery curve, Lewis et al. 2017, Fig S4F
%             (values from the vector drawing, s.e.m. as drawn).
%   D.train   trains of 30 ms steps for 4 s, Lewis et al. 2017, Fig 3B: first and
%             last peak as a fraction of a 300 ms step given 10 s earlier
%             (digitised from the figure, s.e.m. as drawn).
%   D.sust    current left at the end of a 300 ms step, cell attached, -80 mV:
%             0.096 (s.d. 0.040, N = 75), Lewis et al. 2017, Fig S1E.
%===============================================================================

function D = Piezo1_targets()

D.T_rec = 295;                                  % kelvin, room temperature taken as 22 C

%% Fitted
D.decay.V_mV  = [-60 -40 -20 0 20 -60 -40 -20 0 20];
D.decay.tau   = [9.12 11.05 13.95 17.79 24.09 10.75 12.13 14.85 18.64 26.80];   % ms
D.decay.err   = [1.25 1.57 2.00 2.72 3.95 1.17 1.67 2.05 3.20 3.24];            % ms, as plotted
D.decay.set   = [1 1 1 1 1 2 2 2 2 2];          % 1 Wu et al. 2017, 2 Zheng et al. 2019
D.decay.names = {'Wu et al. 2017','Zheng et al. 2019'};

D.deact.tau   = 4.0;                            % ms
D.deact.tol   = 0.2*4.0;                        % ms, 20 %: no error bar is drawn

D.half.delay  = 0.140;                          % second
D.half.tol    = 0.2*0.140;                      % second, 20 %

%% Checked only
D.human.tau   = 8.6;                            % ms
D.human.sem   = 0.4;                            % ms

D.zheng80.tau = 11.9;                           % ms
D.zheng80.sem = 0.6;                            % ms

D.wufar.V_mV  = [-100 -80 40 60 80 100];
D.wufar.tau   = [8.5 8 33 48 70 95];            % ms, read by eye
D.wufar.tol   = 0.2*D.wufar.tau;                % ms, 20 %

D.rec.delay   = [0.01 0.02 0.04 0.08 0.16 0.32 0.64 0.64 1.28 2.56 5.12 10.24 20.48 40.96];   % second
D.rec.value   = [0.1679 0.3163 0.4213 0.4393 0.5142 0.5052 0.5532 0.5982 0.6327 0.6627 0.7616 0.8261 0.8861 0.9745];
D.rec.sem     = [0.0495 0.0570 0.0630 0.0615 0.0585 0.0600 0.0555 0.0765 0.0645 0.0690 0.0750 0.0600 0.0540 0.0630];

D.train.f_Hz  = [0.5 1 2 5 10];
D.train.first     = [0.7804 0.7473 0.8579 0.8895 0.9345];
D.train.first_sem = [0.062 0.073 0.083 0.083 0.106];
D.train.last      = [0.6351 0.4981 0.4002 0.2943 0.2717];
D.train.last_sem  = [0.075 0.070 0.062 0.042 0.055];

D.sust.value  = 0.096;
D.sust.sd     = 0.040;
end
