%===============================================================================
% Action potential duration at 90 % repolarisation
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
%
%   a = apd90_ms(t, Vm_mV, u)
%   t      : time, second
%   Vm_mV  : membrane potential, millivolt
%   u      : upstroke time of the beat, second (upstroke_times.m)
%   a      : APD90 in ms, from the upstroke through -20 mV to 90 % repolarisation
%            within 1.9 s; the reference level is the minimum in the 0.3 s before
%            the upstroke. NaN if the beat does not repolarise in the window.
%===============================================================================

function a = apd90_ms(t, V, u)

a  = nan;
w  = t >= u & t <= u+1.9;
tt = t(w);
VV = V(w);
[Vp,ip] = max(VV);
Vmin    = min(V(t >= u-0.3 & t <= u));
if isempty(Vmin)
    Vmin = min(V);
end
V90     = Vp - 0.9*(Vp-Vmin);
jr      = ip-1+find(VV(ip:end-1) >= V90 & VV(ip+1:end) < V90,1);
if ~isempty(jr)
    a = 1e3*(tt(jr) + (V90-VV(jr))*(tt(jr+1)-tt(jr))/(VV(jr+1)-VV(jr)) - u);
end
end
