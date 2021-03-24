function [vegetation, obu_dif] = ObuFunc(vegetation, p, ic, il, obu_val)
% DESCRIPTION:
% This is the function to solve for the Obukhov length (obu). For the current
% estimate of the Obukhov length (obu_val), calculate u* , T*, and q* and then
% the new length (obu). The function value is the change in Obukhov length and
% equals zero when the the Obukhov length does not change value between iterations.

% ARGUMENTS:

%     p               % Patch index for CLM g/l/c/p hierarchy
%     ic              % Aboveground layer index
%     il              % Sunlit (1) or shaded (2) leaf index
%     obu_val         % Input value for Obukhov length (m)

% mlcanopy_inst

% LOCAL VARIABLES:
%     c1                          % Parameter to calculate beta_neutral
%     beta_neutral                % Neutral value for beta = u*/u(h), ~0.3-0.5
%     beta                        % Value of u*/u(h) at canopy top
%     dt                          % Height below canopy top (dt = ztop - zdisp)
%     zeta                        % Monin-Obukhov stability parameter (z-d)/L
%     psim                        % psi function for momentum
%     psic                        % psi function for scalars
%     zlog                        % log((zref-zdisp)/(ztop-zdisp))
%     tvstar                      % Virtual potential temperature scale (K)
%     obu_dif                     % Difference in Obuhkov length (m)
%---------------------------------------------------------------------

%Input
% % turb_type = 2;
% ztop    = vegetation.mlcanopyinst.ztop;           % Canopy height (m)
% z0mg    = vegetation.mlcanopyinst.z0mg;           % Roughness length of ground (m)
% lai     = vegetation.canopy.lai;            % Leaf area index of canopy (m2/m2)
% sai     = vegetation.mlcanopyinst.sai;            % Stem area index of canopy (m2/m2)
% Lc      = vegetation.mlcanopyinst.Lc;             % Canopy density length scale (m)
% zref    = vegetation.mlcanopyinst.zref;           % Reference height (m)
% uref    = vegetation.mlcanopyinst.uref;           % Wind speed at reference height (m/s)
% rhomol  = vegetation.mlcanopyinst.rhomol;         % Molar density at reference height (mol/m3)
% qref    = vegetation.mlcanopyinst.qref;           % Specific humidity at reference height (kg/kg)
% qaf      = vegetation.mlcanopyinst.qaf;           % Specific humidity at canopy top (kg/kg)
% thref   = vegetation.mlcanopyinst.thref;          % Atmospheric potential temperature (K)
% thvref  = vegetation.mlcanopyinst.thvref;         % Atmospheric virtual potential temperature (K)
% taf     = vegetation.mlcanopyinst.taf;            % Air temperature at canopy top (K)
% 
% %Output
% ustar   = vegetation.mlcanopyinst.ustar;          % Friction velocity (m/s)
% tstar   = vegetation.mlcanopyinst.tstar;          % Temperature scale (K)
% qstar   = vegetation.mlcanopyinst.qstar;          % Water vapor scale (kg/kg)
% gah     = vegetation.mlcanopyinst.gah;            % Aerodynamic conductance for a scalar above canopy (mol/m2/s)
% obu     = obu;            % Obukhov length (m)
% % obu_gah = vegetation.mlcanopyinst.obu_gah;        % Obukhov length used for gah (m)
% PrSc    = vegetation.mlcanopyinst.PrSc;           % Prandtl (Schmidt) number at canopy top
% zdisp   = vegetation.mlcanopyinst.zdisp;          % Displacement height (m)
% uaf     = vegetation.mlcanopyinst.uaf;            % Wind speed at canopy top (m/s)

% p           = vegetation.mlcanopyinst.p;      % Patch index for CLM g/l/c/p hierarchy
% vegetation.pftcon.beta_neutral_max = beta_neutral_max;
% vegetation.pftcon.cr = cr; 

% Use this current value of Obukhov length
vegetation.mlcanopyinst.obu(p) = obu_val;

% Prevent near-zero value of Obukhov length

if (abs(vegetation.mlcanopyinst.obu(p)) <= 0.1) 
    vegetation.mlcanopyinst.obu(p) = 0.1;
end
    
    % Determine neutral value of beta
%     
% % switch(turb_type)
% %    case (1 & 3)
c1 = (vegetation.physcon.vkc / log((vegetation.mlcanopyinst.ztop(p) + vegetation.mlcanopyinst.z0mg(p))/vegetation.mlcanopyinst.z0mg(p)))^2;
beta_neutral = min(sqrt(c1 + vegetation.pftcon.cr*(vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p))), vegetation.pftcon.beta_neutral_max);
% %    case (2)
% %         beta_neutral = vkc / 2;
% % end

% Calculate beta = u* / u(h) for current Obukhov length

% call GetBeta (beta_neutral, Lc(p)/obu(p), beta)

LcL = vegetation.mlcanopyinst.Lc(p)/vegetation.mlcanopyinst.obu(p);
[beta] = GetBeta (beta_neutral, LcL);


% Displacement height: dt is height below canopy top (dt = ztop - zdisp)

dt = beta^2.* vegetation.mlcanopyinst.Lc(p);
dt = dt.* (1  - exp(-0.25.*(vegetation.canopy.lai(p)+vegetation.mlcanopyinst.sai(p))/beta^2));
dt = min(vegetation.mlcanopyinst.ztop(p), dt);
vegetation.mlcanopyinst.zdisp(p) = vegetation.mlcanopyinst.ztop(p) - dt;

if ((vegetation.mlcanopyinst.zref(p) - vegetation.mlcanopyinst.zdisp(p)) < 0.)
    error('ERROR: CanopyTurbulenceMod: ObuFunc error, zdisp height > zref');
end

% Calculate turbulent Prandlt number (Pr) and Schmidt number (Sc) at canopy
% top for current Obukhov length. Only need Pr because Sc = Pr.

% % switch(turb_type)
% %     case (1&3)

[vegetation.mlcanopyinst.PrSc(p)] = GetPrSc (beta_neutral, vegetation.pftcon.beta_neutral_max, vegetation.mlcanopyinst.Lc(p)/vegetation.mlcanopyinst.obu(p));
% %     case (2)
% %         PrSc(p) = 1;
% % end

% Calculate the stability functions psi for momentum and scalars. The
% returned function values (psim, psic) are the Monin-Obukhov psi functions
% and additionally include the roughness sublayer psi_hat functions.
% These are evaluated at the reference height and at the canopy height.
% Here, limit Obukhov length based on values of zeta so that extreme
% cases are excluded.

zeta = (vegetation.mlcanopyinst.zref(p) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
if (zeta >= 0 )
    zeta = min(1 , max(zeta,0.01));
else
    zeta = max(-2 , min(zeta,-0.01)); %
end
vegetation.mlcanopyinst.obu(p) = (vegetation.mlcanopyinst.zref(p) - vegetation.mlcanopyinst.zdisp(p)) / zeta;

% call GetPsiRSL (zref(p), ztop(p), dt, beta, zeta, obu(p), PrSc(p), psim, psic);
[psim, psic] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zref(p), vegetation.mlcanopyinst.ztop, dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));

% Friction velocity

zlog = log((vegetation.mlcanopyinst.zref(p)-vegetation.mlcanopyinst.zdisp(p)) / (vegetation.mlcanopyinst.ztop(p)-vegetation.mlcanopyinst.zdisp(p)));
vegetation.mlcanopyinst.ustar(p) = vegetation.mlcanopyinst.uref(p).* vegetation.physcon.vkc / (zlog + psim);

% Wind speed at canopy height

vegetation.mlcanopyinst.uaf(p) = vegetation.mlcanopyinst.ustar(p) / beta;

% Temperature scale

% % % % both added here from CanopyTurbulence.m
% % % vegetation.mlcanopyinst.taf(p) = vegetation.mlcanopyinst.tair(p,vegetation.canopy.ntop(p));
% % % vegetation.mlcanopyinst.qaf(p) = vegetation.physcon.mmh2o / vegetation.physcon.mmdry.* vegetation.mlcanopyinst.eaf(p) / (vegetation.mlcanopyinst.pref(p) - (1  - vegetation.physcon.mmh2o/vegetation.physcon.mmdry).* vegetation.mlcanopyinst.eaf(p));

vegetation.mlcanopyinst.tstar(p) = (vegetation.mlcanopyinst.thref(p) - vegetation.mlcanopyinst.taf(p)).* vegetation.physcon.vkc / (zlog + psic);

% Water vapor scale - use units of specific humidity (kg/kg)

vegetation.mlcanopyinst.qstar(p) = (vegetation.mlcanopyinst.qref(p) - vegetation.mlcanopyinst.qaf(p)).* vegetation.physcon.vkc / (zlog + psic);

% Aerodynamic conductance to canopy height

vegetation.mlcanopyinst.gah(p) = vegetation.mlcanopyinst.rhomol(p).* vegetation.physcon.vkc .* vegetation.mlcanopyinst.ustar(p) / (zlog + psic);
vegetation.mlcanopyinst.obu_gah(p) = vegetation.mlcanopyinst.obu(p);

% Calculate new Obukhov length (m)

tvstar = vegetation.mlcanopyinst.tstar(p) + 0.61 .* vegetation.mlcanopyinst.thref(p).* vegetation.mlcanopyinst.qstar(p);
vegetation.mlcanopyinst.obu(p) = vegetation.mlcanopyinst.ustar(p)^2.* vegetation.mlcanopyinst.thvref(p) / (vegetation.physcon.vkc .* vegetation.physcon.grav .* tvstar);

% Change in Obukhov length (m)

obu_dif = vegetation.mlcanopyinst.obu(p) - obu_val(p); %obu_val;

end