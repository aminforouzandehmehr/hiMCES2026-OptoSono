%===============================================================================
% Times at which the membrane potential rises through -20 mV
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
%
%   u = upstroke_times(t, Vm_mV)
%   t      : time, second
%   Vm_mV  : membrane potential, millivolt
%   u      : upstroke times, second (linear interpolation between samples)
%===============================================================================

function u = upstroke_times(t, V)

j = find(V(1:end-1) < -20 & V(2:end) >= -20);
u = t(j) + (-20-V(j)).*(t(j+1)-t(j))./(V(j+1)-V(j));
end
