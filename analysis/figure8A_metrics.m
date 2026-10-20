%===============================================================================
% Measurements of the Fig. 8A protocol
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
%
%   M = figure8A_metrics(light, stim, nostim, off, stimOn, stimDur, fid)
%
%   light   : run with the three light pulses at 0, 2 and 4 s (0 to 4.9 s)
%   stim    : continuation with the mechanical pulse (4.9 to 9.5 s)
%   nostim  : continuation without the mechanical pulse
%   off     : continuation of the Piezo1-off cell with the mechanical pulse
%   stimOn  : onset of the mechanical pulse, second
%   stimDur : duration of the mechanical pulse, second
%   fid     : optional file identifier; the report is also written there
%
% A light pulse captures if an upstroke through -20 mV follows within 0.1 s of
% its onset. The mechanical pulse captures if an upstroke falls between its
% onset and 0.1 s after its end and the run without it has none there.
%===============================================================================

function M = figure8A_metrics(light, stim, nostim, off, stimOn, stimDur, fid)

if nargin < 7
    fid = [];
end
win = [stimOn stimOn+stimDur+0.1];             % second

%% Light pulses
ul   = upstroke_times(light.t,light.Vm_mV);
pulses = [0 2 4];                              % second
M.light_captured = 0;
for p = pulses
    M.light_captured = M.light_captured + any(ul >= p & ul <= p+0.1);
end
M.peak_IChR2_AF   = min(light.I_ChR2);
M.APD90_light_ms  = apd90_ms(light.t,light.Vm_mV,ul(2));   % the light beat at 2 s

%% Mechanical pulse and controls
us = upstroke_times(stim.t,stim.Vm_mV);
un = upstroke_times(nostim.t,nostim.Vm_mV);
uo = upstroke_times(off.t,off.Vm_mV);
js = find(us >= win(1) & us <= win(2),1);
M.mech_captured        = ~isempty(js) && ~any(un >= win(1) & un <= win(2));
M.nostim_in_window     = any(un >= win(1) & un <= win(2));
M.Piezo1off_in_window  = any(uo >= win(1) & uo <= win(2));
M.nostim_upstrokes_s   = un(:)';
w = stim.t >= stimOn & stim.t <= stimOn+stimDur;
M.peak_IPiezo1_in_AF   = min(stim.I_Piezo1(w));
M.peak_IPiezo1_out_AF  = max(stim.I_Piezo1(w));
if isempty(js)
    M.latency_ms       = nan;
    M.APD90_mech_ms    = nan;
    M.peak_Cai_mech_nM = nan;
else
    M.latency_ms       = 1e3*(us(js)-stimOn);
    M.APD90_mech_ms    = apd90_ms(stim.t,stim.Vm_mV,us(js));
    M.peak_Cai_mech_nM = max(stim.Cai_nM(stim.t >= stimOn & stim.t < stimOn+1.5));
end
M.peak_Cai_light_nM    = max(light.Cai_nM(light.t >= 4));

%% Report
say(fid,'Light pulses (10 ms, 0.15 mW/mm2, 0.5 Hz): %d of 3 captured; peak I_ChR2 %.2f A/F\n',M.light_captured,M.peak_IChR2_AF);
say(fid,'Mechanical pulse at %.1f s for %.1f s: captured %d, latency %.1f ms\n',stimOn,stimDur,M.mech_captured,M.latency_ms);
say(fid,'  peak I_Piezo1 during the pulse: %.2f A/F inward, %+.2f A/F outward\n',M.peak_IPiezo1_in_AF,M.peak_IPiezo1_out_AF);
say(fid,'Controls: no mechanical pulse, upstroke in the window %d (next upstrokes at %s s);\n',M.nostim_in_window,mat2str(round(M.nostim_upstrokes_s,3)));
say(fid,'          Piezo1 off with the mechanical pulse, upstroke in the window %d\n',M.Piezo1off_in_window);
say(fid,'APD90: light beat at 2 s %.0f ms, mechanical beat %.0f ms\n',M.APD90_light_ms,M.APD90_mech_ms);
say(fid,'Peak Ca2+: light beat at 4 s %.0f nM, mechanical beat %.0f nM\n',M.peak_Cai_light_nM,M.peak_Cai_mech_nM);
end

function say(fid, varargin)
fprintf(varargin{:});
if ~isempty(fid)
    fprintf(fid,varargin{:});
end
end
