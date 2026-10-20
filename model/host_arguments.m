%===============================================================================
% The arguments of the published hiMCES2026 that follow time and Y, at their
% standard values: active contraction, no drug, no ischemia, no isoproterenol.
% The names are those of the published MasterCompute_hiMCES.m.
%
%   A = host_arguments(stimFlag, freq, classicVSoptic)
%   stimFlag         : 0 no electrical stimulus, 1 electrical stimulus
%   freq             : pacing frequency in hertz
%   classicVSoptic   : 0 standard solution (Na 151, K 5.4, Ca 1.8 mM), 37 C
%                      1 Tyrode solution   (Na 135, K 5.4, Ca 1.33 mM), 21 C
%                      X.T_override changes the temperature and leaves the bath.
%
% Use as
%   [dY, dati] = hiMCES2026_ChR2_Piezo1(time, Y, A{:}, X);
%   [dY, dati] = hiMCES2026(time, Y(1:33), A{:});   % upstream/hiMCES2026.m
%===============================================================================

function A = host_arguments(stimFlag, freq, classicVSoptic)

contrflag       = 1;                 % 0 isometric, 1 active contraction, 2 no contraction
GAC             = ones(1,25);        % conductance scales; ones for the standard model
tDrugApplication = 1e5;              % second, later than the end of any run here

% Reduction and scaling factors (1 = no effect)
INaFRedMed      = 1;
INaLRedMed      = 1;
INaKRedMed      = 1;
ICaLRedMed      = 1;
IKrRedMed       = 1;
IKsRedMed       = 1;
INaCaRedMed     = 1;
IfRedMed        = 1;
SERCARedMed     = 1;
SERCAKmRedMed   = 1;
RyRRedMed       = 1;
RyRPoRedMed     = 1;
ItoRedMed       = 1;

ischemiaFlag    = 0;                 % 0 none, 1 SEV1, 2 SEV2
myFlags         = [0 0 0 0 0 0 0];   % [HyperK, Acidosis, IKatp, INaK, IpCa, Jup, CE]
IKatpSelectionString = 'Kazbanov';   % 'Ledezma' or 'Kazbanov'
fnp             = 0;                 % physioxia
z               = 1;                 % physioxia
flv             = 1;
fmyo            = ones(1,15);
gma             = 22;                % standard value
ISOFlag         = 0;                 % 0 OFF, 1 ON

A = {contrflag, GAC, stimFlag, freq, classicVSoptic, tDrugApplication, ...
     INaFRedMed, INaLRedMed, INaKRedMed, ICaLRedMed, IKrRedMed, IKsRedMed, INaCaRedMed, ...
     IfRedMed, SERCARedMed, SERCAKmRedMed, RyRRedMed, RyRPoRedMed, ItoRedMed, ...
     ischemiaFlag, myFlags, IKatpSelectionString, fnp, z, flv, fmyo, gma, ISOFlag};
end
