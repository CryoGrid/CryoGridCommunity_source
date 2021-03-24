function [vegetation] = CanopyTurbulence (vegetation, p, niter) %(p, niter, soilstate_inst, temperature_inst, energyflux_inst, waterflux_inst, vegetation)

% DESCRIPTION:
% Canopy turbulence, aeorodynamic conductances, and wind/temperature/water vapor
% profiles using above- and within-canopy coupling with a roughness sublayer
% (RSL) parameterization

% USES:
%     turb_type %clm_varctl
%     use_scalars %clm_varctl
%     mmh2o %clm_varcon
%     mmdry %clm_varcon
%     vegetation.physcon.vkc  %clm_varcon

%     use MathToolsMod, only : hybrid
%     use PatchType, only : patch
%     use SoilStateType, only : soilstate_type
%     use TemperatureType, only : temperature_type
%     use EnergyFluxType, only : energyflux_type
%     use WaterFluxType, only : waterflux_type
%     use CanopyFluxesMultilayerType, only : mlcanopy_type

%     pftcon % pftconMod

% ARGUMENTS:
%     p       % Patch index for CLM g/l/c/p hierarchy
%     niter   % Iteration index
%     soilstate_inst
%     temperature_inst
%     energyflux_inst
%     waterflux_inst
%     mlcanopy_inst

% LOCAL VARIABLES:
%     c                  % Column index for CLM g/l/c/p hierarchy
%     ic                 % Canopy layer index
%     il                 % Sunlit (1) or shaded (2) leaf index
%     iterate            % CLM: 1 = iterate for Obukhov length. 0 = use specified Obukhov length
%     dummy              % Dummy argument
%     lm                 % Length scale for momentum (m)
%     beta               % Value of u*/u(h) at canopy top
%     res                % Resistance (s/m)
%     zu, zl             % Upper and lower heights for within canopy resistances (m)
%     obu0, obu1         % Initial estimates for Obukhov length (m)
%     z0hg               % Roughness length of ground for sensible heat (m)
%     zlog_m, zlog_c     % Log-profiles at ground
%     ustar_g            % Friction velocity at ground (m/s)
%     eta                % Used to limit value for beta/lm
%     beta_over_lm       % beta/lm (1/m)
%     tol                % Accuracy tolerance for Obukhov length (m)
%     zeta               % (z - d)/L
%     phic               % Flux-gradient relationship at (hc - d)/L
%     kc_at_hc           % Scalar eddy diffusivity at canopy top (m2/s)

%%% NEW
%     dt                 % Height below canopy top (dt = vegetation.mlcanopyinst.ztop - zdisp)
%     psim               % psi function for momentum
%     psic               % psi function for scalars
%     psim1, psim2
%     psic1, psic2
%     ga_add
%     sumres

%     zdisp_most, z0m_most, zeta_most, psim_most_z0m, psim_most, zlog_most, ustar_most
%     z0c_most, sh_most, tstar_most, psic_most_z0c, psic_most, ts_most
%---------------------------------------------------------------------

% Input
% p           = vegetation.mlcanopyinst.p;      % Patch index for CLM g/l/c/p hierarchy
% ic          = vegetation.mlcanopyinst.ic;     % Canopy layer index
% il          = vegetation.mlcanopyinst.il;     % Sunlit (1) or shaded (2) leaf index


% vegetation.mlcanopyinst.vegetation.mlcanopyinst.ztop     = vegetation.mlcanopyinst.vegetation.mlcanopyinst.ztop;                % Canopy height (m)
% vegetation.mlcanopyinst.ncan     = vegetation.mlcanopyinst.ncan;                % Number of aboveground layers
% vegetation.canopy.ntop     = vegetation.mlcanopyinst.vegetation.canopy.ntop;                % Index for top leaf layer
% vegetation.canopy.lai      = vegetation.mlcanopyinst.lai;                 % Leaf area index of canopy (m2/m2)
% vegetation.mlcanopyinst.sai      = vegetation.mlcanopyinst.sai;                 % Stem area index of canopy (m2/m2)
% vegetation.mlcanopyinst.zs       = vegetation.mlcanopyinst.zs;                  % Canopy height for scalar concentration and source (m)
% vegetation.mlcanopyinst.co2ref   = vegetation.mlcanopyinst.co2ref;              % Atmospheric CO2 at reference height (umol/mol)
% vegetation.mlcanopyinst.pref     = vegetation.mlcanopyinst.pref;                % Air pressure at reference height (Pa)
% vegetation.mlcanopyinst.rhomol   = vegetation.mlcanopyinst.rhomol;              % Molar density at reference height (mol/m3)
% vegetation.mlcanopyinst.eg       = vegetation.mlcanopyinst.eg;                  % Soil surface vapor pressure (Pa)
% vegetation.mlcanopyinst.tg       = vegetation.mlcanopyinst.tg;                  % Soil surface temperature (K)
% vegetation.mlcanopyinst.zref     = vegetation.mlcanopyinst.zref;                % Reference height (m)
% vegetation.mlcanopyinst.uref     = vegetation.mlcanopyinst.uref;                % Wind speed at reference height (m/s)
% vegetation.mlcanopyinst.thref    = vegetation.mlcanopyinst.thref;               % Atmospheric potential temperature (K)
% vegetation.mlcanopyinst.cpair    = vegetation.mlcanopyinst.cpair;               % Specific heat of air at constant pressure, at reference height (J/mol/K)
%
% vegetation.pftcon.z0mr     = vegetation.pftcon.z0mr;                       % Ratio of momentum roughness length to canopy top height
% vegetation.pftcon.displar  = vegetation.pftcon.displar;                    % Ratio of displacement height to canopy top height
%
% % Output
% vegetation.mlcanopyinst.z0mg     = vegetation.mlcanopyinst.z0mg;                % Roughness length of ground (m)
% vegetation.mlcanopyinst.Lc       = vegetation.mlcanopyinst.Lc;                  % Canopy density length scale (m)
% vegetation.mlcanopyinst.obu      = vegetation.mlcanopyinst.obu;                 % Obukhov length (m)
% vegetation.mlcanopyinst.obu_gah  = vegetation.mlcanopyinst.vegetation.mlcanopyinst.obu_gah;          % Obukhov length used for gah (m)
% vegetation.mlcanopyinst.obuold   = vegetation.mlcanopyinst.obuold;                % Obukhov length from previous iteration
% vegetation.mlcanopyinst.nmozsgn  = vegetation.mlcanopyinst.vegetation.mlcanopyinst.nmozsgn;          % Number of times stability changes sign during iteration
% vegetation.mlcanopyinst.ustar    = vegetation.mlcanopyinst.ustar;                 % Friction velocity (m/s)
% % ????????????? vegetation.mlcanopyinst.tstar    = vegetation.mlcanopyinst.tstar;                 % Temperature scale (K)
% % ????????????? vegetation.mlcanopyinst.qstar    = vegetation.mlcanopyinst.qstar;                 % Water vapor scale (kg/kg)
% vegetation.mlcanopyinst.PrSc     = vegetation.mlcanopyinst.PrSc;                % Prandtl (Schmidt) number at canopy top
% vegetation.mlcanopyinst.zdisp    = vegetation.mlcanopyinst.zdisp;                 % Displacement height (m)
% vegetation.mlcanopyinst.uaf      = vegetation.mlcanopyinst.uaf;                 % Wind speed at canopy top (m/s)
% vegetation.mlcanopyinst.taf      = vegetation.mlcanopyinst.taf;                 % Air temperature at canopy top (K)
% vegetation.mlcanopyinst.eaf      = vegetation.mlcanopyinst.eaf;                 % Vapor pressure at canopy top (Pa)
% vegetation.mlcanopyinst.qaf      = vegetation.mlcanopyinst.qaf;                 % Specific humidity at canopy top (kg/kg)
% vegetation.mlcanopyinst.wind     = vegetation.mlcanopyinst.wind;                % Wind speed profile (m/s)
% vegetation.mlcanopyinst.wind_most= vegetation.mlcanopyinst.wind_most;        % Wind speed profile from MOST (m/s)
% vegetation.mlcanopyinst.tair     = vegetation.mlcanopyinst.tair;                % Air temperature profile (K)
% vegetation.mlcanopyinst.tair_most= vegetation.mlcanopyinst.tair_most;        % Air temperature profile from MOST (K)
% vegetation.mlcanopyinst.eair     = vegetation.mlcanopyinst.eair;                % Vapor pressure profile (Pa)
% vegetation.mlcanopyinst.cair     = vegetation.mlcanopyinst.cair;                % Atmospheric CO2 profile (umol/mol)
% vegetation.mlcanopyinst.gah      = vegetation.mlcanopyinst.gah;                 % Aerodynamic conductance for a scalar above canopy (mol/m2/s)
% vegetation.mlcanopyinst.ga_prof  = vegetation.mlcanopyinst.ga_prof;        % Canopy layer aerodynamic conductance for scalars (mol/m2/s)
% vegetation.mlcanopyinst.cd = vegetation.mlcanopyinst.cd;
% turb_type = vegetation.physcon.turb_type;


% Initialization for first iteration
obuold(p) = vegetation.mlcanopyinst.obu(p);
if (niter == 1)
    obuold(p) = 0;
    vegetation.mlcanopyinst.nmozsgn(p) = 0;
end

% Set roughness length of ground

vegetation.mlcanopyinst.z0mg(p) = 0.01;

% Canopy density length scale

vegetation.mlcanopyinst.Lc(p) = vegetation.mlcanopyinst.ztop(p) / (vegetation.mlcanopyinst.cd .* (vegetation.canopy.lai(p) + vegetation.mlcanopyinst.sai(p)));

% These are not used, but are needed to pass into the hybrid root solver

ic = 0;
il = 0;

% Calculate Obukhov length (obu)

% % switch(turb_type)
% %
% %     case(1&3)

% Calculate Obukhov length (obu) using the function ObuFunc to iterate until
% the change in obu is < tol. do not use the returned value (dummy), and
% instead use the value of obu calculated in the final call to ObuFunc.

obu0 = 100;           % Initial estimate for Obukhov length (m)
obu1 = -100;          % Initial estimate for Obukhov length (m)
tol = 0.1;             % Accuracy tolerance for Obukhov length (m)
func_name = 'ObuFunc';

[vegetation, ~] = hybrid_root(func_name, vegetation, p, ic, il, obu0, obu1, tol); %('CanopyTurbulence', p, ic, il, mlcanopy_inst, ObuFunc, obu0, obu1, tol)

vegetation.mlcanopyinst.obu(p) = vegetation.mlcanopyinst.obu_gah(p);

% % case (4)

% Use CLM calculation of Obukhov length (obu)
%end
% iterate = 1;
% % [vegetation] = ObuFuncCLM (vegetation);
% % vegetation.mlcanopyinst.obu(p) = vegetation.mlcanopyinst.obu_gah(p);
% % end

% Check to see if Obukhov length is changing signs between iterations.
% If too many changes in sign, set it to a near-neutral value.

if (obuold(p)*vegetation.mlcanopyinst.obu(p) < 0 )
    vegetation.mlcanopyinst.nmozsgn(p) = vegetation.mlcanopyinst.nmozsgn(p) + 1;
    obuold(p) = vegetation.mlcanopyinst.obu(p);
end
if (vegetation.mlcanopyinst.nmozsgn(p) >= 4)
    vegetation.mlcanopyinst.obu(p) = -1000;
    % %         switch(turb_type)
    % %             case(1:3)
    vegetation = ObuFunc (vegetation, p, ic, il, vegetation.mlcanopyinst.obu(p)); %ML: why dummy?
    vegetation.mlcanopyinst.obu(p) = vegetation.mlcanopyinst.obu_gah(p);
    % %             case (4)
    % % %                 iterate = 0;
    % % %                 [vegetation] = ObuFuncCLM (vegetation);
    % % %                 vegetation.mlcanopyinst.obu(p) = vegetation.mlcanopyinst.obu_gah(p);
    % %         end
end

% Above-canopy wind profile: wind speed is defined at zs

% % switch(turb_type)
% %     case(1&2&3)
dt = vegetation.mlcanopyinst.ztop(p) - vegetation.mlcanopyinst.zdisp(p);
beta = vegetation.mlcanopyinst.ustar(p) / vegetation.mlcanopyinst.uaf(p);


for ic = vegetation.canopy.ntop(p)+1:vegetation.mlcanopyinst.ncan(p)
    zeta = (vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
    [psim, ~] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zs(p,ic), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu, vegetation.mlcanopyinst.PrSc);
    zlog_m = log((vegetation.mlcanopyinst.zs(p,ic)-vegetation.mlcanopyinst.zdisp(p)) / (vegetation.mlcanopyinst.ztop(p)-vegetation.mlcanopyinst.zdisp(p)));
    vegetation.mlcanopyinst.wind(p,ic) = vegetation.mlcanopyinst.ustar(p) / vegetation.physcon.vkc.* (zlog_m + psim);
end

% % end

% Above-canopy aerodynamic conductances: these are defined between
% vegetation.mlcanopyinst.zs(i) and vegetation.mlcanopyinst.zs(i+1). Note the special case for the top layer to the
% reference height.

% % switch(turb_type)
% %     case(1&2&3)

for ic = vegetation.canopy.ntop(p)+1:vegetation.mlcanopyinst.ncan(p)-1
    zeta = (vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
    [~, psic1] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zs(p,ic), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
    zeta = (vegetation.mlcanopyinst.zs(p,ic+1) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
    [~, psic2] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zs(p,ic+1), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
    % equivalent to: -psi_c_z2 + psi_c_z1 + psi_c_rsl_z2 - psi_c_rsl_z1
    psic = psic2 - psic1;
    zlog_c = log((vegetation.mlcanopyinst.zs(p,ic+1)-vegetation.mlcanopyinst.zdisp(p)) / (vegetation.mlcanopyinst.zs(p,ic)-vegetation.mlcanopyinst.zdisp(p)));
    vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.rhomol(p).* vegetation.physcon.vkc.* vegetation.mlcanopyinst.ustar(p) / (zlog_c + psic);
end

ic = vegetation.mlcanopyinst.ncan(p);
zeta = (vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
[~, psic1] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zs(p,ic), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
zeta = (vegetation.mlcanopyinst.zref(p) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
[~, psic2] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zref(p), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
psic = psic2 - psic1;
zlog_c = log((vegetation.mlcanopyinst.zref(p)-vegetation.mlcanopyinst.zdisp(p)) / (vegetation.mlcanopyinst.zs(p,ic)-vegetation.mlcanopyinst.zdisp(p)));
vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.rhomol(p).* vegetation.physcon.vkc.* vegetation.mlcanopyinst.ustar(p) / (zlog_c + psic);
% % end

%Within-canopy wind profile: wind speed is defined at vegetation.mlcanopyinst.zs

beta = vegetation.mlcanopyinst.ustar(p) / vegetation.mlcanopyinst.uaf(p);
lm = 2. .* beta^3. .* vegetation.mlcanopyinst.Lc(p);

% % switch(turb_type)
% %     case (1&2)
%eta = min (beta/lm*vegetation.mlcanopyinst.ztop(p), 10 );
eta = min (beta/lm*vegetation.mlcanopyinst.ztop(p), 20);
% %     case (3&4)
% %         eta = 3 ;
% % end

% --------------------------------------------------
% ---------------------------------------------------
% -----------------------------------------------------

beta_over_lm = eta / vegetation.mlcanopyinst.ztop(p);
%   beta_over_lm = beta / lm

for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p) %Fortran ic = 1:ncan
    vegetation.mlcanopyinst.wind(p,ic) = vegetation.mlcanopyinst.uaf(p).* exp((vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.ztop(p)).* beta_over_lm);
    vegetation.mlcanopyinst.wind(p,ic) = max(vegetation.mlcanopyinst.wind(p,ic), 0.1);
end

% Scalar eddy diffusivity at canopy top

% % switch(turb_type)
% %
% %     case (1&2&3)

kc_at_hc = 2 .* beta^3.* vegetation.mlcanopyinst.Lc(p).* vegetation.mlcanopyinst.ustar(p) / vegetation.mlcanopyinst.PrSc(p);

% %     case (4)
% %
% %         zeta = (vegetation.mlcanopyinst.ztop(p) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
% %         if (zeta < zetah)               % very unstable (zeta < zetah)
% %             phic = 0.9.* vegetation.physcon.vkc^1.333 .* (-zeta)^(-0.333 );
% %         elseif (zeta < 0 )          % unstable (zetah <= zeta < 0)
% %             phic = 1  / sqrt(1  - 16 .* zeta);
% %         elseif (zeta <=  1 )        % stable (0 <= zeta <= 1)
% %             phic = 1  + 5 .* zeta;
% %         else                               % very stable (zeta > 1)
% %             phic = 5  + zeta;
% %         end
% %         kc_at_hc = vegetation.physcon.vkc.* vegetation.mlcanopyinst.ustar(p).* (vegetation.mlcanopyinst.ztop(p) - vegetation.mlcanopyinst.zdisp(p)) / phic;
% %
% % end

% Within-canopy aerodynamic conductances: these are defined between
% vegetation.mlcanopyinst.zs(i) and vegetation.mlcanopyinst.zs(i+1). Note the special case for the top layer to the
% canopy height.

for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)-1 %Fortran ic = 1:ncan
    zl = vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.ztop(p);
    zu = vegetation.mlcanopyinst.zs(p,ic+1) - vegetation.mlcanopyinst.ztop(p);
    %      res = vegetation.mlcanopyinst.PrSc(p) / (beta.* vegetation.mlcanopyinst.ustar(p)).* (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
    res = 1  / (kc_at_hc.* beta_over_lm).* (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm));
    vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.rhomol(p) / res;
end

% Special case for top layer: conductance to top of canopy

ic = vegetation.canopy.ntop(p);
zl = vegetation.mlcanopyinst.zs(p,ic) - vegetation.mlcanopyinst.ztop(p);
zu = vegetation.mlcanopyinst.ztop(p) - vegetation.mlcanopyinst.ztop(p);
%   res = vegetation.mlcanopyinst.PrSc(p) / (beta.* vegetation.mlcanopyinst.ustar(p)).* (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm))
res = 1.  / (kc_at_hc.* beta_over_lm).* (exp(-zl*beta_over_lm) - exp(-zu*beta_over_lm));
vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.rhomol(p) / res;

% Now add additional resistance to first atmospheric layer

if (vegetation.mlcanopyinst.ncan(p) == vegetation.canopy.ntop(p))
    ga_add = vegetation.mlcanopyinst.gah(p);
elseif (vegetation.mlcanopyinst.ncan(p) > vegetation.canopy.ntop(p))
    % %     switch(turb_type)
    % %         case (1&2&3)
    zeta1 = (vegetation.mlcanopyinst.ztop(p) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
    %             call GetPsiRSL (vegetation.mlcanopyinst.ztop(p), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p), psim1, psic1);
    [~, psic1] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.ztop(p), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta1, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
    zeta2 = (vegetation.mlcanopyinst.zs(p,ic+1) - vegetation.mlcanopyinst.zdisp(p)) / vegetation.mlcanopyinst.obu(p);
    %             call GetPsiRSL (vegetation.mlcanopyinst.zs(p,ic+1), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p), psim2, psic2);
    [~, psic2] = GetPsiRSL(vegetation, vegetation.mlcanopyinst.zs(p,ic+1), vegetation.mlcanopyinst.ztop(p), dt, beta, zeta2, vegetation.mlcanopyinst.obu(p), vegetation.mlcanopyinst.PrSc(p));
    psic = psic2 - psic1;
    zlog_c = log((vegetation.mlcanopyinst.zs(p,ic+1)-vegetation.mlcanopyinst.zdisp(p)) / (vegetation.mlcanopyinst.ztop(p)-vegetation.mlcanopyinst.zdisp(p)));
    ga_add = vegetation.mlcanopyinst.rhomol(p).* vegetation.physcon.vkc.* vegetation.mlcanopyinst.ustar(p) / (zlog_c + psic);
    % %     end
end

vegetation.mlcanopyinst.ga_prof(p,ic) = 1  / (1  / vegetation.mlcanopyinst.ga_prof(p,ic) + 1  / ga_add);

% Make sure above-canopy aerodynamic resistances sum to 1/vegetation.mlcanopyinst.gah
sumres = 0;
for ic = vegetation.canopy.ntop(p)+1:vegetation.mlcanopyinst.ncan(p)
    sumres = sumres + 1. / vegetation.mlcanopyinst.ga_prof(p,ic);
end
sumres = sumres + 1. / ga_add;

if (abs(1./sumres - vegetation.mlcanopyinst.gah(p)) > 1.e-06 )
    disp(1./sumres - vegetation.mlcanopyinst.gah(p));
    disp(' ERROR: CanopyTurbulenceMod: above-canopy aerodynamic conductance error');
end

% Resistance at ground

%   zl = 0  - vegetation.mlcanopyinst.ztop(p)
%   zu = vegetation.mlcanopyinst.zs(p,1) - vegetation.mlcanopyinst.ztop(p)
%   res = vegetation.mlcanopyinst.PrSc(p) / (beta.* vegetation.mlcanopyinst.ustar(p)).* (exp(-zl*beta/lm) - exp(-zu*beta/lm))
%   res = vegetation.mlcanopyinst.ga_prof(p,1)

%SEBAS: BONAN commented out these lines above and replaced it with the
%stuff below - might be worth considering to change back? This could be
%passed on to the GORUND module?

%this is taken care of in GROUND module now!
z0hg = 0.1 .* vegetation.mlcanopyinst.z0mg(p);
zlog_m = log(vegetation.mlcanopyinst.zs(p,2)/vegetation.mlcanopyinst.z0mg(p));
zlog_c = log(vegetation.mlcanopyinst.zs(p,2)/z0hg);
ustar_g = vegetation.mlcanopyinst.wind(p,2).* vegetation.physcon.vkc / zlog_m;
ustar_g = max(ustar_g, 0.01);
res = zlog_c / (vegetation.physcon.vkc.* ustar_g);
vegetation.mlcanopyinst.ga_prof(p,1) = vegetation.mlcanopyinst.rhomol(p) / res; %was p,0

if (zlog_m < 0 || zlog_c < 0 )
    disp(' ERROR: CanopyTurbulenceMod: soil roughness error');
end



% Convert resistance (s/m) to conductance (mol/m2/s)
% Limit resistances to < 500 s/m

for ic = 1:vegetation.mlcanopyinst.ncan(p)   %ic = 0:vegetation.mlcanopyinst.ncan(p)
    res = min(vegetation.mlcanopyinst.rhomol(p)/vegetation.mlcanopyinst.ga_prof(p,ic), 500.);
    vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.rhomol(p) / res;
end

% Values at ground
vegetation.mlcanopyinst.wind(p,1) = 0.; %p,0
%CHANGED SEBAS
%vegetation.mlcanopyinst.tair(p,1) = vegetation.mlcanopyinst.tg(p); %p,0
vegetation.mlcanopyinst.tair(p,1) = vegetation.mlcanopyinst.tair(p,2);
%END CHANGE SEBAS!
%vegetation.mlcanopyinst.eair(p,1) = vegetation.mlcanopyinst.eg(p); %p,0

% Calculate within-canopy scalar profiles for temperature and vapor pressure
% using an implicit solution

% call ScalarProfile (p, niter, soilstate_inst, temperature_inst, energyflux_inst, waterflux_inst, mlcanopy_inst);
[vegetation] = scalar_profile (vegetation, p);

% Commented out in FORTRAN original
% % Calculate within-canopy scalar profiles for temperature, vapor pressure, and CO2
% % using an iterative coupling
% %
% % Temperature: scalar is T*cp (J/mol) and flux is J/m2/s
% % Water vapor: scalar is e/P (mol/mol) and flux is mol/m2/s
% % CO2: scalar is c (umol/mol) and flux is umol/m2/s
%
% [cval] = ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), tref(p), vegetation.mlcanopyinst.ga_prof(p,:), vegetation.mlcanopyinst.sh_prof(p,:), 1/vegetation.mlcanopyinst.cpair(p), vegetation.mlcanopyinst.tair_old(p,:), vegetation.mlcanopyinst.tair(p,:));
% % call ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), tref(p), vegetation.mlcanopyinst.ga_prof(p,:), sh_prof(p,:), &
% % 1 /vegetation.mlcanopyinst.cpair(p), tair_old(p,:), vegetation.mlcanopyinst.tair(p,:))
%
% [cval] = ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), eref(p), vegetation.mlcanopyinst.ga_prof(p,:), vegetation.mlcanopyinst.et_prof(p,:), vegetation.mlcanopyinst.pref(p), vegetation.mlcanopyinst.eair_old(p,:), vegetation.mlcanopyinst.eair(p,:));
% % call ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), eref(p), vegetation.mlcanopyinst.ga_prof(p,:), et_prof(p,:), &
% % vegetation.mlcanopyinst.pref(p), eair_old(p,:), vegetation.mlcanopyinst.eair(p,:))
%
% [cval] = ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), vegetation.mlcanopyinst.co2ref(p), vegetation.mlcanopyinst.ga_prof(p,:), vegetation.mlcanopyinst.fc_prof(p,:),1 , vegetation.mlcanopyinst.cair_old(p,:), vegetation.mlcanopyinst.cair(p,:));
% % call ScalarProfileIterative (vegetation.mlcanopyinst.ncan(p), vegetation.mlcanopyinst.co2ref(p), vegetation.mlcanopyinst.ga_prof(p,:), fc_prof(p,:), &
% % 1 , vegetation.mlcanopyinst.cair_old(p,:), cair(p,:))

for ic = 1:vegetation.mlcanopyinst.ncan(p)      %ic = 0:vegetation.mlcanopyinst.ncan(p)
    vegetation.mlcanopyinst.cair(p,ic) = vegetation.mlcanopyinst.co2ref(p);
end

% if (~ use_scalars)
% for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p) %Fortran ic = 1:ncan
%     vegetation.mlcanopyinst.wind(p,ic) = vegetation.mlcanopyinst.wind(p,vegetation.canopy.ntop(p));
%     vegetation.mlcanopyinst.tair(p,ic) = vegetation.mlcanopyinst.tair(p,vegetation.canopy.ntop(p));
%     vegetation.mlcanopyinst.eair(p,ic) = vegetation.mlcanopyinst.eair(p,vegetation.canopy.ntop(p));
% end
% end

vegetation.mlcanopyinst.taf(p) = vegetation.mlcanopyinst.tair(p,vegetation.canopy.ntop(p));
vegetation.mlcanopyinst.eaf(p) = vegetation.mlcanopyinst.eair(p,vegetation.canopy.ntop(p));
vegetation.mlcanopyinst.qaf(p) = vegetation.physcon.mmh2o / vegetation.physcon.mmdry.* vegetation.mlcanopyinst.eaf(p) / (vegetation.mlcanopyinst.pref(p) - (1  - vegetation.physcon.mmh2o/vegetation.physcon.mmdry).* vegetation.mlcanopyinst.eaf(p));

% vegetation.mlcanopyinst.wind and temperature profiles using MOST

% vegetation.mlcanopyinst.wind_most(p,0) = -999;
% vegetation.mlcanopyinst.tair_most(p,0) = -999;
vegetation.mlcanopyinst.wind_most(p,1) = -999.; %p,0
vegetation.mlcanopyinst.tair_most(p,1) = -999.; %p,0

zdisp_most = vegetation.pftcon.displar(p).* vegetation.mlcanopyinst.ztop(p);
z0m_most = vegetation.pftcon.z0mr(p).* vegetation.mlcanopyinst.ztop(p);
z0c_most = z0m_most;

zeta_most = (vegetation.mlcanopyinst.zref(p) - zdisp_most) / vegetation.mlcanopyinst.obu(p);
psim_most_z0m = psi_m_monin_obukhov(z0m_most/vegetation.mlcanopyinst.obu(p));

psim_most = -psi_m_monin_obukhov(zeta_most) + psim_most_z0m;

zlog_most = log((vegetation.mlcanopyinst.zref(p)-zdisp_most) / z0m_most);
ustar_most = vegetation.mlcanopyinst.uref(p).* vegetation.physcon.vkc / (zlog_most + psim_most);

sh_most = -vegetation.mlcanopyinst.cpair(p).* (vegetation.mlcanopyinst.tair(p,vegetation.canopy.ntop(p)+1) - vegetation.mlcanopyinst.tair(p,vegetation.canopy.ntop(p))).* vegetation.mlcanopyinst.ga_prof(p,vegetation.canopy.ntop(p));
tstar_most = -sh_most / (vegetation.mlcanopyinst.rhomol(p).* vegetation.mlcanopyinst.cpair(p).* ustar_most);
psic_most_z0c = psi_c_monin_obukhov(z0c_most/vegetation.mlcanopyinst.obu(p));

psic_most = -psi_c_monin_obukhov(zeta_most) + psic_most_z0c;

zlog_most = log((vegetation.mlcanopyinst.zref(p)-zdisp_most) / z0c_most);
ts_most = vegetation.mlcanopyinst.thref(p) - tstar_most / vegetation.physcon.vkc.* (zlog_most + psic_most);

for ic = vegetation.mlcanopyinst.ncan(p):-1:2 %fortran 1
    if (vegetation.mlcanopyinst.zs(p,ic) > zdisp_most)
        zeta_most = (vegetation.mlcanopyinst.zs(p,ic) - zdisp_most) / vegetation.mlcanopyinst.obu(p);
        psim_most = -(psi_m_monin_obukhov(zeta_most)) + psim_most_z0m;
        psic_most = -(psi_c_monin_obukhov(zeta_most)) + psic_most_z0c;
        zlog_most = log((vegetation.mlcanopyinst.zs(p,ic)-zdisp_most) / z0m_most);
        vegetation.mlcanopyinst.wind_most(p,ic) = ustar_most / vegetation.physcon.vkc.* (zlog_most + psim_most);
        zlog_most = log((vegetation.mlcanopyinst.zs(p,ic)-zdisp_most) / z0c_most);
        vegetation.mlcanopyinst.tair_most(p,ic) = ts_most + tstar_most / vegetation.physcon.vkc.* (zlog_most + psic_most);
    else
        vegetation.mlcanopyinst.wind_most(p,ic) = -999.;
        vegetation.mlcanopyinst.tair_most(p,ic) = -999.;
    end
end

end



