%===============================================================================
% Integrates the two gates of the modified Piezo1 model at a fixed voltage through a
% stimulus given as straight segments.
%
%   S = channel_run(X, V_mV, T, tb, Fb, tout)
%
%   X      : options from ChR2_Piezo1_opts
%   V_mV   : membrane potential in millivolt, constant
%   T      : temperature in kelvin
%   tb, Fb : break points of the stimulus: times in seconds and the tension that
%            the activation gate senses there, in kPa. Straight lines between.
%   tout   : times in seconds at which the result is wanted
%   S      : t, ma, h and g = ma*h, the open fraction, at those times
%
% The run starts from rest: ma at its steady state for Fb(1), and h = 1.
% The integration is cut at every break point of the stimulus.
%===============================================================================

function S = channel_run(X, V_mV, T, tb, Fb, tout)

V    = 1e-3*V_mV;                                   % volt
opt  = odeset('RelTol',1e-9,'AbsTol',1e-12);
tout = tout(:);
ma   = zeros(size(tout));
h    = zeros(size(tout));

[~,~,r] = Piezo1_channel(V,0,1,Fb(1),T,X);
y       = [r.ma_inf; 1];

for k = 1:numel(tb)-1
    a   = tb(k);
    b   = tb(k+1);
    if b <= a
        continue
    end
    Fa  = Fb(k);
    Fs  = (Fb(k+1)-Fb(k))/(b-a);
    sol = ode45(@(t,y) rhs(t,y,V,T,X,a,Fa,Fs),[a b],y,opt);
    if k < numel(tb)-1
        w = tout >= a & tout < b;
    else
        w = tout >= a & tout <= b;
    end
    if any(w)
        yy    = deval(sol,tout(w));
        ma(w) = yy(1,:)';
        h(w)  = yy(2,:)';
    end
    y = sol.y(:,end);
end

S.t  = tout;
S.ma = ma;
S.h  = h;
S.g  = ma.*h;
end

function dy = rhs(t, y, V, T, X, a, Fa, Fs)
F        = Fa + Fs*(t-a);                           % kPa
[dma,dh] = Piezo1_channel(V,y(1),y(2),F,T,X);
dy       = [dma; dh];
end
