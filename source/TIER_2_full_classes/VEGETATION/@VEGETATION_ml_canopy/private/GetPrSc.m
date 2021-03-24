function [Pr] = GetPrSc (beta_neutral, beta_neutral_max, LcL)
%
% DESCRIPTION:
% Calculate turbulent Prandlt number (Pr) and Schmidt number (Sc) at canopy
% top for current Obukhov length
%
% ARGUMENTS:
% Input:
%     beta_neutral_max = 0.35 ;
%     beta_neutral = 0.3;
%     LcL = 1;

% Output:
%Pr  %turbulent Prandlt number (Pr)
%Sc                 % Schmidt number at canopy top

% LOCAL VARIABLES:
Prn = 0.5;        % Neutral value for Pr
Prvr = 0.3;       % Magnitude of variation of Pr with stability
Prsc = 2.0;       % Scale of variation of Pr with stability

Scn = Prn;          % Neutral value for Sc
Scvr = Prvr;        % Magnitude of variation of Sc with stability
Scsc = Prsc;        % Scale of variation of Sc with stability
%---------------------------------------------------------------------

Pr = Prn + Prvr.* tanh(Prsc*LcL);
Pr = (1  - beta_neutral/beta_neutral_max).* 1  + beta_neutral/beta_neutral_max*Pr;

Sc = Scn + Scvr.* tanh(Scsc*LcL);
Sc = (1  - beta_neutral/beta_neutral_max).* 1  + beta_neutral/beta_neutral_max*Sc;
end