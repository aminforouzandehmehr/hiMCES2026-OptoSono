%===============================================================================
% Version 1.0.0 (software DOI: 10.5281/zenodo.23430167)
% Author: Amin Forouzandehmehr
%
%   ok = verify_release()          equations, Fig. 8A, Fig. S5 (about 1 min)
%   ok = verify_release('full')    also conditioning of both cells and the
%                                  full channel checks (about 3 min more)
%
% Writes results/verification_report.txt and returns true if every check
% passes. Reference values are in reference_outputs/ and initial_conditions/.
%
% Tolerances for the Fig. 8A traces are those of the development checks:
% Vm 0.05 mV, Cai 0.5 nM, I_Piezo1 2e-5 A/F, I_ChR2 0.02 A/F. States and
% channel measurements: relative 1e-6. On the computer that produced the
% references the differences are expected to be zero.
%===============================================================================

function ok = verify_release(level)

if nargin == 0
    level = 'quick';
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'model'));
addpath(fullfile(here,'analysis'));
addpath(fullfile(here,'calibration'));
addpath(fullfile(here,'upstream'));
outdir = fullfile(here,'results');
if ~isfolder(outdir)
    mkdir(outdir);
end
fid = fopen(fullfile(outdir,'verification_report.txt'),'w');
cleanup = onCleanup(@() fclose(fid));
log = @(varargin) both(fid,varargin{:});
log('verify_release(''%s''), MATLAB %s, %s\n\n',level,version,datestr(now,'yyyy-mm-dd HH:MM'));
ok = true;

%% 1) Equations
log('1) EQUATIONS\n');
ok = check_equations(here,log) && ok;

%% 2) Settings of Fig. 8A
log('\n2) SETTINGS\n');
ok = check_settings(here,log) && ok;

%% 3) Fig. 8A and controls
log('\n3) FIGURE 8A AND CONTROLS\n');
ok = check_figure8A(here,log) && ok;

%% 4) Fig. S5 and channel measurements
log('\n4) FIGURE S5 AND CHANNEL MEASUREMENTS\n');
ok = check_figureS5(here,log) && ok;

%% 5) Conditioning and full channel checks
if strcmp(level,'full')
    log('\n5) CONDITIONING\n');
    ok = check_conditioning(here,log) && ok;
    log('\n6) FULL CHANNEL CHECKS (calibration/check_Piezo1.m)\n');
    M = check_Piezo1();
    R = load(fullfile(here,'reference_outputs','channel_measurements.mat'));
    ok = report(log,ok,relmax(M,R.M) <= 1e-6, ...
        sprintf('check_Piezo1 measurements = reference, largest relative difference %.1e; report in results/check_report.txt',relmax(M,R.M)));
end

log('\n%s\n',verdict(ok));
end

%===============================================================================
% Checks
%===============================================================================

function ok = check_equations(here, log)
ok = true;
F  = 96485.3415;      % coulomb_per_mole
Vc = 8800;            % micrometre_cube
Cm = 9.87109e-11;     % farad

% 400 sampled states: host states with 5 % scatter, Vm from -90 to +50 mV,
% random gates and channel fractions; the random stream is fixed
S  = load(fullfile(here,'verification','start_states.mat'));
rs = RandStream('mt19937ar','Seed',7);
n  = 400;
Ys = zeros(39,n);
ts = 2*rand(rs,1,n);
for i = 1:n
    y       = S.Y33(:,mod(i-1,size(S.Y33,2))+1).*(1 + 0.05*(2*rand(rs,33,1)-1));
    y(1)    = 1e-3*(-90 + 140*rand(rs));
    c       = rand(rs,3,1);
    c       = c/sum(c)*rand(rs);
    Ys(:,i) = [y; rand(rs); 2*rand(rs); c; rand(rs)];
end

% a) channels off: the published host in every derivative and output
X0 = ChR2_Piezo1_opts('gp_scale',0);
mx = 0;
mo = 0;
for i = 1:n
    for cvo = 0:1
        A         = host_arguments(mod(i,2),1,cvo);
        y         = Ys(:,i);
        y(36:39)  = 0;
        [d39,o39] = hiMCES2026_ChR2_Piezo1(ts(i),y,A{:},X0);
        [d33,o33] = hiMCES2026(ts(i),y(1:33),A{:});
        mx        = max(mx,max(abs(d39(1:33)-d33(:))));
        mo        = max(mo,max(abs(o39(1:48)-o33(:)')));
    end
end
ok = report(log,ok,mx == 0 && mo == 0,sprintf( ...
    'a) Piezo1 and ChR2 off = upstream/hiMCES2026.m: largest difference %g in dY, %g in outputs (%d states, both baths)',mx,mo,n));

% b) gating option 'released' = Piezo1.m version 1.1
Tset = {[],[],304};
cvo  = [0 1 1];
Tk   = [310 294 304];
e    = 0;
for i = 1:n
    k     = mod(i-1,3)+1;
    A     = host_arguments(mod(i,2),1,cvo(k));
    X     = ChR2_Piezo1_opts('gating','released','accounting','published','T_override',Tset{k});
    y     = Ys(:,i);
    [d,o] = hiMCES2026_ChR2_Piezo1(ts(i),y,A{:},X);
    [dma,dXn,I,J] = piezo1_v11(y(1),y(34),y(35),o(21),Tk(k));
    e     = max([e, abs(d(34)-dma)/max(abs(dma),1e-12), abs(d(35)-dXn)/max(abs(dXn),1e-12), ...
                 abs(o(49)-I)/max(abs(I),1e-15), abs(o(50)-J)/max(abs(J),1e-15)]);
end
[same,miss,nst] = same_lines(fullfile(here,'upstream','Piezo1.m'),fullfile(here,'model','hiMCES2026_ChR2_Piezo1.m'));
ok = report(log,ok,e < 1e-12 && same,sprintf( ...
    'b) option gating = released: block = transcription of Piezo1.m v1.1 (relative %.1e); %d of %d statements of upstream/Piezo1.m occur unchanged, the other 2 differ only by the added options',e,nst-miss,nst));

% c) modified gating: coupled block = model/Piezo1_channel.m (used by the calibration)
e = 0;
for i = 1:n
    k     = mod(i-1,3)+1;
    A     = host_arguments(mod(i,2),1,cvo(k));
    X     = ChR2_Piezo1_opts('T_override',Tset{k});
    y     = Ys(:,i);
    [d,o] = hiMCES2026_ChR2_Piezo1(ts(i),y,A{:},X);
    [dma,dh,r] = Piezo1_channel(y(1),y(34),y(35),o(51),Tk(k),X);
    e     = max([e, abs(d(34)-dma)/max(abs(dma),1e-12), abs(d(35)-dh)/max(abs(dh),1e-12), ...
                 abs(o(52)-r.ma_inf), abs(o(54)-r.tau_inact)/r.tau_inact, abs(o(60)-r.tau_rec)/r.tau_rec]);
end
ok = report(log,ok,e < 1e-12,sprintf( ...
    'c) modified gating: coupled model = isolated channel (Piezo1_channel.m) at 310, 294 and 304 K, relative %.1e',e));

% d) split ion accounting
X   = ChR2_Piezo1_opts('T_override',304,'touchF',0.056,'touchN',1e9,'touchDur',10,'touchFreq',0.01);
Xn  = X;
Xn.gp_scale = 0;
Buf_C  = 0.25-(70e-3);
Kbuf_C = 0.001;
s1  = 0;
s2  = 0;
s3  = 0;
s4  = 0;
for i = 1:n
    A     = host_arguments(0,0.5,1);
    y     = Ys(:,i);
    [d,o] = hiMCES2026_ChR2_Piezo1(ts(i),y,A{:},X);
    dn    = hiMCES2026_ChR2_Piezo1(ts(i),y,A{:},Xn);
    s1    = max(s1,abs((dn(1)-d(1)) - o(49))/max(abs(o(49)),1e-9));
    bufc  = 1.0/(1.0 + Buf_C*Kbuf_C/(y(3)+Kbuf_C)^2);
    if o(49) < 0
        caIn = bufc*X.f_Ca*o(49)*Cm/(2.0*Vc*F*1e-18);
        naIn = (1.0-X.f_Ca)*o(49)*Cm/(F*Vc*1e-18);
        s2   = max(s2,abs((d(3)-dn(3)) + caIn)/max(abs(caIn),1e-12));
        s3   = max(s3,abs((d(18)-dn(18)) + naIn)/max(abs(naIn),1e-12));
    else
        s4   = max(s4,max(abs(d([3 18]) - dn([3 18]))));
    end
end
ok = report(log,ok,s1 < 1e-9 && s2 < 1e-4 && s3 < 1e-9 && s4 == 0,sprintf( ...
    'd) accounting = split: Vm carries I_Piezo1 (%.1e); inward current shared Ca2+ (%.1e) and Na+ (%.1e); outward current in no balance (%g)',s1,s2,s3,s4));

% e) round voltages: undefined only at exactly 0 V, where the published host is undefined
X      = ChR2_Piezo1_opts('T_override',304);
bad    = 0;
badAt0 = 0;
nanpub = 0;
finAdd = true;
y      = Ys(:,1);
y(34:39) = [0.5; 0.5; 0.1; 0.1; 0.2; 0.5];
for c = 0:1
    A = host_arguments(0,1,c);
    for v = -1000:600
        y(1)   = v*1e-4;
        [d,o]  = hiMCES2026_ChR2_Piezo1(0.5,y,A{:},X);
        isBad  = ~(all(isfinite(d)) && all(isfinite(o)));
        bad    = bad + isBad;
        badAt0 = badAt0 + (isBad && v == 0);
        finAdd = finAdd && all(isfinite(d(34:39))) && all(isfinite(o(49:58)));
    end
    y(1)   = 0;
    d33    = hiMCES2026(0.5,y(1:33),A{:});
    nanpub = nanpub + any(~isfinite(d33));
end
ok = report(log,ok,bad == badAt0 && badAt0 == nanpub && finAdd,sprintf( ...
    'e) -100 to +60 mV in 0.1 mV steps: undefined values %d, all at 0 V as in the published host (%d); added states finite',bad,nanpub));

% f) options: conductance scale, applied stimulus
A  = host_arguments(0,0.5,1);
e6 = 0;
e7 = 0;
for i = 1:40
    y      = Ys(:,i);
    [~,o1] = hiMCES2026_ChR2_Piezo1(0.5,y,A{:},ChR2_Piezo1_opts('T_override',304));
    [~,o5] = hiMCES2026_ChR2_Piezo1(0.5,y,A{:},ChR2_Piezo1_opts('T_override',304,'gp_scale',5));
    e6     = max(e6,abs(o5(49)-5*o1(49))/max(abs(o1(49)),1e-15));
    Xt     = ChR2_Piezo1_opts('T_override',304,'touchF',0.03,'touchT0',1,'touchN',1);
    [~,a]  = hiMCES2026_ChR2_Piezo1(0.98,y,A{:},Xt);
    [~,b]  = hiMCES2026_ChR2_Piezo1(1.02,y,A{:},Xt);
    [~,c]  = hiMCES2026_ChR2_Piezo1(1.06,y,A{:},Xt);
    e7     = max([e7, abs(a(51)-a(21)), abs(b(51)-b(21)-0.03), abs(c(51)-c(21))]);
end
ok = report(log,ok,e6 < 1e-12 && e7 < 1e-15,sprintf( ...
    'f) I_Piezo1 scales with gp_scale (%.1e); the applied stimulus adds to the sensed stress only while on (%.1e kPa)',e6,e7));

% g) native stress filter
tau = 0.15;
Xs  = ChR2_Piezo1_opts('T_override',304,'tau_trans',tau);
Xt  = ChR2_Piezo1_opts('T_override',304,'tau_trans',tau,'touchF',0.03,'touchT0',1,'touchN',1);
X0  = ChR2_Piezo1_opts('T_override',304);
e1  = 0;
e2  = 0;
e3  = 0;
e4  = 0;
for i = 1:40
    y40     = [Ys(:,i); 0.2*mod(i*0.6180339887498949,1)];
    [d0,o0] = hiMCES2026_ChR2_Piezo1(0.5,Ys(:,i),A{:},X0);
    [ds,os] = hiMCES2026_ChR2_Piezo1(0.5,y40,A{:},Xs);
    [~,ot]  = hiMCES2026_ChR2_Piezo1(1.02,y40,A{:},Xt);
    e1      = max(e1,abs(ds(40) - (os(21)-y40(40))/tau));
    e2      = max(e2,abs(os(51) - y40(40)));
    e3      = max(e3,abs(ot(51) - y40(40) - 0.03));
    de      = hiMCES2026_ChR2_Piezo1(0.5,[Ys(:,i); o0(21)],A{:},Xs);
    e4      = max(e4,max(abs(de(1:39)-d0(1:39))./max(abs(d0(1:39)),1e-12)));
end
ok = report(log,ok,e1 < 1e-15 && e2 < 1e-15 && e3 < 1e-15 && e4 < 1e-9,sprintf( ...
    'g) state 40 = first order filter of the native stress (%.1e); gate senses it (%.1e) plus the unfiltered stimulus (%.1e); other states unchanged (%.1e)',e1,e2,e3,e4));
end

function ok = check_settings(here, log)
ok = true;
c  = {'Piezo1','Piezo1off'};
for i = 1:2
    L  = load(fullfile(here,'initial_conditions',sprintf('OptoSono_31C_%s.mat',c{i})));
    X  = OptoSono_settings(c{i});
    Xs = L.X;
    X.optT0  = Xs.optT0;
    X.optN   = Xs.optN;
    ok = report(log,ok,isequaln(X,Xs),sprintf( ...
        '%s: OptoSono_settings = options stored with the conditioned state (pulse timing excluded); %d states',c{i},numel(L.Y)));
end
end

function ok = check_figure8A(here, log)
ok   = true;
cvo  = 1;
dt   = 2e-4;
IC1  = load(fullfile(here,'initial_conditions','OptoSono_31C_Piezo1.mat'));
IC0  = load(fullfile(here,'initial_conditions','OptoSono_31C_Piezo1off.mat'));
X1   = OptoSono_settings('Piezo1');
X0   = OptoSono_settings('Piezo1off');
XL1  = X1;
XL1.optN = 3;
XL0  = X0;
XL0.optN = 3;
mech = @(X,F) setfield(setfield(setfield(X,'touchF',F),'touchT0',5.0),'touchN',1); %#ok<SFLD>

tic
L1 = run_segments(XL1,IC1.Y,0,4.9,dt,cvo);
S1 = run_segments(mech(X1,0.336),L1.Yend,4.9,9.5,dt,cvo);
N1 = run_segments(X1,L1.Yend,4.9,9.5,dt,cvo);
L0 = run_segments(XL0,IC0.Y,0,4.9,dt,cvo);
N0 = run_segments(X0,L0.Yend,4.9,9.5,dt,cvo);
S0 = run_segments(mech(X0,1.344),L0.Yend,4.9,9.5,dt,cvo);
log('   six runs in %.0f s\n',toc);

% traces against reference_outputs/Figure8A_traces.csv (every 1 ms)
T   = readtable(fullfile(here,'reference_outputs','Figure8A_traces.csv'));
k1  = 1:5:numel(L1.t)-1;
k2  = 1:5:numel(S1.t);
cc  = @(a,b,f) [a.(f)(k1); b.(f)(k2)];
tol = struct('Vm_mV',0.05,'Cai_nM',0.5,'I_ChR2',0.02,'I_Piezo1',2e-5);
sets = {'Piezo1_stim',L1,S1; 'Piezo1_nostim',L1,N1; 'Piezo1off_nostim',L0,N0};
cols = {'Vm_mV','Vm_mV'; 'Cai_nM','Cai_nM'; 'I_ChR2','IChR2_AF'; 'I_Piezo1','IPiezo1_AF'};
for s = 1:3
    for c = 1:4
        name = [sets{s,1} '_' cols{c,2}];
        if ~ismember(name,T.Properties.VariableNames)
            continue
        end
        dmax = max(abs(cc(sets{s,2},sets{s,3},cols{c,1}) - T.(name)));
        ok   = report(log,ok,dmax <= tol.(cols{c,1}),sprintf('%-28s max |difference| %.3g (tolerance %g)',name,dmax,tol.(cols{c,1})));
    end
end

% end states against reference_outputs/Figure8A_end_states.mat
E  = load(fullfile(here,'reference_outputs','Figure8A_end_states.mat'));
pr = {'light_Piezo1',L1; 'stim_Piezo1',S1; 'nostim_Piezo1',N1; 'light_Piezo1off',L0; 'nostim_Piezo1off',N0};
for i = 1:size(pr,1)
    r  = relvec(pr{i,2}.Yend,E.(pr{i,1}));
    ok = report(log,ok,r <= 1e-6,sprintf('end state %-17s largest relative difference %.1e',pr{i,1},r));
end

% events and controls
M  = figure8A_metrics(L1,S1,N1,S0,5.0,0.5);
ok = report(log,ok,M.light_captured == 3,'light: 3 of 3 pulses captured');
ok = report(log,ok,M.mech_captured,sprintf('mechanical pulse 0.336 kPa captured, latency %.1f ms',M.latency_ms));
ok = report(log,ok,~M.nostim_in_window,'control without the mechanical pulse: no upstroke from 5.0 to 5.6 s');
ok = report(log,ok,~M.Piezo1off_in_window,'Piezo1-off control with a 1.344 kPa pulse: no upstroke from 5.0 to 5.6 s');

% threshold bracket
for F = [0.336*2^(5/2)/16 0.168]
    Q  = run_segments(mech(X1,F),L1.Yend,4.9,9.5,dt,cvo);
    u  = upstroke_times(Q.t,Q.Vm_mV);
    hit = any(u >= 5.0 & u <= 5.6);
    ok = report(log,ok,hit == (F > 0.15),sprintf('threshold scan: %.3f kPa, beat %d', F, hit));
end
end

function ok = check_figureS5(here, log)
ok = true;
X  = ChR2_Piezo1_opts();
D  = Piezo1_targets();
Dc = D;
Dc.decay.V_mV = -100:5:100;
Mc = measure_Piezo1(X,Dc,false);
C  = readtable(fullfile(here,'reference_outputs','FigureS5_curve.csv'));
r  = max(abs(Mc.decay_tau(:)-C.decay_tau_ms)./C.decay_tau_ms);
ok = report(log,ok,r <= 1e-6,sprintf('Fig. S5 curve, 41 voltages: largest relative difference %.1e',r));
M  = measure_Piezo1(X,D,true);
R  = load(fullfile(here,'reference_outputs','channel_measurements.mat'));
r  = relmax(M,R.M);
ok = report(log,ok,r <= 1e-6,sprintf('fitted targets and checks: largest relative difference %.1e',r));
log('   reported values (failed checks retained): human %.2f ms, sustained %.4f, +40 to +100 mV %s ms, recovery at 1.28 s %.3f\n', ...
    M.human_tau,M.sust,mat2str(round(M.wufar_tau(3:6),1)),M.rec(9));
end

function ok = check_conditioning(here, log)
ok = true;
c  = {'Piezo1','Piezo1off'};
for i = 1:2
    st = condition_OptoSono(c{i});
    L  = load(fullfile(here,'initial_conditions',sprintf('OptoSono_31C_%s.mat',c{i})));
    r  = relvec(st.Y,L.Y);
    ok = report(log,ok,st.ok && st.T == L.conditioning.duration_s && r <= 1e-6,sprintf( ...
        '%s: converged %d after %g s (reference %g s); state vs initial_conditions, largest relative difference %.1e', ...
        c{i},st.ok,st.T,L.conditioning.duration_s,r));
end
end

%===============================================================================
% Local functions
%===============================================================================

function [dma, dXn, I, J] = piezo1_v11(V, ma, Xn, frc, T)
% Transcription of Piezo1.m version 1.1 (upstream/Piezo1.m).
d      = (frc/0.056) * 5;
mpi    = 1.0/(1+exp(-(d-2.5)/0.25));
dma    = (mpi-ma)/0.1;
q10    = 2^((310-T)/10);
Xs_inf = 2.0/(1.0+exp((-V*1000.0-0.0)/16.0));
al_Xs  = 1/sqrt(1.0+exp((-60.0+V*1000.0)/50.0));
be_Xs  = 900.0/(1.0+exp((+171.0-V*1000.0)/72.0));
ta_Xs  = 1e-3 * q10 * 30.0*al_Xs*be_Xs/100.0;
dXn    = (Xs_inf-Xn)/ta_Xs;
I      = 0.3039*ma*Xn*(V-0);
J      = I * 0.5*1.1626e-01;
end

function [same, miss, nst] = same_lines(released, model)
% Every statement of the released file must occur in the model, comments and
% spaces removed. Two statements carry the added options (conductance scale,
% stress input); they must equal the released ones once the option is removed.
strip  = @(L) regexprep(regexprep(L,'%.*$',''),'\s','');
R      = strip(regexp(fileread(released),'\r?\n','split'));
M      = strip(regexp(fileread(model),'\r?\n','split'));
R      = R(~cellfun(@isempty,R));
Mplain = strrep(strrep(strrep(M,'*X.gp_scale',''),'frc_Pie','frc'),'X.F_ref','0.056');
asIs   = ismember(R,M);
plain  = ismember(R,Mplain);
nst    = numel(R);
miss   = sum(~asIs);
same   = all(asIs | plain) && miss == 2 && all(startsWith(R(~asIs),{'g_p=','d='}));
end

function r = relvec(a, b)
r = max(abs(a(:)-b(:))./max(abs(b(:)),1e-12));
end

function r = relmax(A, B)
% Largest relative difference over all numeric fields of two structures.
r  = 0;
fn = fieldnames(B);
for i = 1:numel(fn)
    a = A.(fn{i});
    b = B.(fn{i});
    if isnumeric(b)
        r = max(r,max(abs(a(:)-b(:))./max(abs(b(:)),1e-12)));
    end
end
end

function ok = report(log, ok, pass, msg)
log('%s  %s\n',verdict(pass),msg);
ok = ok && pass;
end

function s = verdict(pass)
if pass
    s = 'PASS';
else
    s = 'FAIL';
end
end

function both(fid, varargin)
fprintf(varargin{:});
fprintf(fid,varargin{:});
end
