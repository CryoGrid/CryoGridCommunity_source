function [vegetation] = scalar_profile (vegetation, p)

% Compute scalar profiles for temperature and water vapor using an implicit
% solution. The boundary condition is the above-canopy scalar values at the
% reference height. The vegetation and ground temperature and fluxes are
% calculated as part of the implicit solution.

%SEBAS change: the ground (+ground surface) temperature as assumed to be fixed, but the
%fluxes (E, H) are calculated instead - these are then passed on to the
%GROUND module to update the energy content of the first cell

% -----------------------------------------------------------------------------------
% Input
%   vegetation.physcon.dtime_sub                  ! Model time step (s)
%   p                   ! Index for grid point to process
%   vegetation.vegetation.physcon.tfrz        ! Freezing point of water (K)
%   vegetation.physcon.hvap        ! Latent heat of evaporation (J/kg)
%   vegetation.physcon.mmh2o       ! Molecular mass of water (kg/mol)
%   vegetation.mlcanopyinst.rhomol      ! Molar density (mol/m3)
%   vegetation.mlcanopyinst.cpair       ! Specific heat of air at constant pressure (J/mol/K)
%   vegetation.mlcanopyinst.pref        ! Atmospheric pressure (Pa)
%   vegetation.mlcanopyinst.thref       ! Potential temperature at reference height (K)
%   vegetation.mlcanopyinst.qref        ! Water vapor at reference height (mol/mol)
%   vegetation.mlcanopyinst.        ! Index for top level
%   surfvar.ntop        ! Index for top leaf layer
%   vegetation.canopy.nbot        ! Index for bottom leaf layer
%   vegetation.soilvar.nsoi        ! First canopy layer is soil
%   vegetation.mlcanopyinst.zw          ! Canopy height at layer interfaces (m)
%   vegetation.canopy.dpai        ! Layer plant area index (m2/m2)
%   surfvar.fwet        ! Fraction of plant area index that is wet
%   surfvar.fdry        ! Fraction of plant area index that is green and dry
%   vegetation.flux.fracsun     ! Sunlit fraction of canopy layer
%   vegetation.flux.fracsha     ! Shaded fraction of canopy layer
%   vegetation.mlcanopyinst.nleaf       ! Number of leaf types (sunlit and shaded)

%   vegetation.mlcanopyinst.gbh         ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
%   leafvar.gbv         ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
%   leafvar.gs          ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
%   leafvar.cpleaf      ! Leaf heat capacity (J/m2 leaf/K)
%   leafvar.rnleaf      ! Leaf net radiation (W/m2 leaf)
%   vegetation.physcon.thk          ! Soil thermal conductivity (W/m/K)
%   vegetation.soilvar.dz          ! Soil layer depth (m)
%   vegetation.soilvar.tsoi        ! Soil temperature (K)
%   soilvar.resis       ! Soil evaporative resistance (s/m)
%   vegetation.mlcanopyinst.rhg         ! Relative humidity of airspace at soil surface (fraction)
%   vegetation.mlcanopyinst.rnsoi       ! Net radiation at ground (W/m2)
%   vegetation.mlcanopyinst.tair_old    ! Air temperature profile for previous timestep (K)
%   vegetation.mlcanopyinst.eair_old    ! Water vapor profile for previous timestep (mol/mol)
%   vegetation.mlcanopyinst.tveg_old    ! Vegetation temperature profile for previous timestep (K)
%   vegetation.mlcanopyinst.ga_prof     ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
% Output
%   vegetation.mlcanopyinst.tair        ! Air temperature profile (K)
%   vegetation.mlcanopyinst.eair        ! Water vapor profile (mol/mol)
%   vegetation.mlcanopyinst.tveg        ! Vegetation temperature profile (K)
%   vegetation.mlcanopyinst.shsoi       ! Ground sensible heat flux, ground (W/m2)
%   vegetation.mlcanopyinst.etsoi       ! Ground evaporation flux (mol H2O/m2/s)
%   vegetation.mlcanopyinst.gsoi        ! Soil heat flux (W/m2)
%   vegetation.mlcanopyinst.tg          ! Soil surface temperature (K)
%   vegetation.mlcanopyinst.shleaf       ! Leaf sensible heat flux (W/m2 leaf)
%   vegetation.mlcanopyinst.etveg       ! Leaf evapotranspiration flux (mol H2O/m2 leaf/s)
%   vegetation.mlcanopyinst.stleaf       ! Leaf storage heat flux (W/m2 leaf)
%   vegetation.mlcanopyinst.shair       ! Canopy air sensible heat flux (W/m2)
%   vegetation.mlcanopyinst.etair       ! Canopy air water vapor flux (mol H2O/m2/s)
%   vegetation.mlcanopyinst.stair       ! Canopy air storage heat flux (W/m2)
% -----------------------------------------------------------------------------------

% isun = 1;      % Sunlit leaf index
% isha = 2;     % Shaded leaf index
% 
% % Latent heat of vaporization (J/mol)
% 
gleaf_sh = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);            % Leaf conductance for sensible heat (mol/m2/s)
gleaf_et = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);            % Leaf conductance for water vapor (mol/m2/s)
alpha = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);               % Coefficient for leaf temperature
beta = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);               
delta = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);               
dqsat = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);               
qsat_term = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);
qsat = zeros(0);
heatcap = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);
avail_energy = zeros(vegetation.mlcanopyinst.ncan,vegetation.mlcanopyinst.nleaf);

% rho_dz_over_dt = zeros(0);     %Intermediate calculation for canopy air storage

a1 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
b11 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
b12 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
c1 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
d1 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature

a2 = zeros(vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
b21 = zeros(vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
b22 = zeros(vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
c2 = zeros(vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
d2 = zeros(vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
% ga_prof_ic_minus_one = zeros(0);     %Special case for ic-1 = 0: use gs0 not ga_prof

% ainv = zeros(0);
% binv = zeros(0);     %"a" and "b" elements of 2x2 matrix to invert
% cinv = zeros(0);
% dinv = zeros(0);     %"c" and "d" elements of 2x2 matrix to invert
% det = zeros(0);     %Determinant of 2x2 matrix

%(0:nlevcan) ??
e11 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
e12 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
f1 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air temperature
e21 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
e22 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)
f2 = zeros(1,vegetation.mlcanopyinst.ncan);     %Coefficient for canopy air water vapor (mole fraction)

%! These are needed only for conservation checks

storage_sh = zeros(1,vegetation.mlcanopyinst.ncan);     %Heat storage flux in air (W/m2)
storage_et = zeros(1,vegetation.mlcanopyinst.ncan);     %Water vapor storage flux in air (mol/m2/s)
stveg = zeros(1,vegetation.mlcanopyinst.ncan);     %Leaf storage heat flux (W/m2)
shsrc = zeros(1,vegetation.mlcanopyinst.ncan);     %Leaf sensible heat flux (W/m2)
etsrc = zeros(1,vegetation.mlcanopyinst.ncan);     %Leaf water vapor flux (mol/m2/s)
% stveg_leaf = zeros(0);     %Leaf storage heat flux from LeafTemperatureMod (W/m2)
% shsrc_leaf = zeros(0);     %Leaf sensible heat flux from LeafTemperatureMod (W/m2)
% etsrc_leaf = zeros(0);     %Leaf water vapor flux from LeafTemperatureMod (mol/m2/s)
% error = zeros(0);     %Energy imbalance (W/m2)
% sum_src = zeros(0);     %Sum of source flux over leaf layers
% sum_storage = zeros(0);     %Sum of storage flux over leaf layers


c=p;

% % % FORTRAN ScalarProfile l. 1410
% % % Get current step (counter)
% % nstep = get_nstep();

vegetation.params.dtime = vegetation.params.dtime_sub;

lambda = LatVap(vegetation.mlcanopyinst.tref(p), vegetation);

isun = vegetation.params.sun; % Array index for sunlit leaf
isha = vegetation.params.sha; % Array index for shaded leaf


% --------------------------------------------------------------------------
% Terms for ground temperature, which is calculated from the energy balance:
% Rn0 - H0 - lambda*E0 - G = 0
% and is rewritten as:
% T0 = alpha0*T(1) + beta0*q(1) + delta0
% --------------------------------------------------------------------------

% %CHANGED SEBAS - use fixed ground surface temperature = temperature of first GROUND grid
% %cell instead
% %[esat, desat] = Satvap (vegetation.mlcanopyinst.tair_old(p,1)); % Pa
% [esat, desat] = Satvap (vegetation.mlcanopyinst.tg(p)); % Pa
% %This is CLM 4.5, compare Eq. 5.73 in CLM 4.5 manual
% %rhg is fraction of saturation of air in the ground cell 
% vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential(c,1) .* 9.81 ./ ((8.3145e+03 ./ 18.016) .*vegetation.mlcanopyinst.tg(p)));
% %added beta factor from CLM 4.5, this is really what seems to supress E -
% %make dependent on field capacity!!
% betaCLM4_5 = 1 +  double(vegetation.soilvar.h2osoi_vol(c,1)<0.5) .* (-1 +  0.25 .* (1-(cos(pi() .* vegetation.soilvar.h2osoi_vol(c,1) ./ 0.5))).^2);
% %this must be multiplied with gs0 to obtain total conductiance
% %end CHANGE SEBAS
% 
% vegetation.mlcanopyinst.eg(p) = esat .* vegetation.mlcanopyinst.rhg(p);
% 
% % qsat0 = esat / vegetation.mlcanopyinst.pref(p);                               % Pa -> mol/mol
% % dqsat0 = desat / vegetation.mlcanopyinst.pref(p);                             % Pa -> mol/mol
% 
% %soilresis is CLM 5, but this needs to depend on water content as
% %described in the manual - seems to be constant here, so E is not regulated
% %down - replaced by CLM 4.5 formulation, although CLM 5 is more physical
% % gsw = 1. / vegetation.mlcanopyinst.soilresis(c);                       %! Soil conductance for water vapor: s/m -> m/s
% % gsw = gsw * vegetation.mlcanopyinst.rhomol(p);                          %! m/s -> mol H2O/m2/s
% % gs0 = vegetation.mlcanopyinst.ga_prof(p,1) * gsw / (vegetation.mlcanopyinst.ga_prof(p,1) + gsw) ; %p,0  %! Total conductance
% 
% gs0 = vegetation.mlcanopyinst.ga_prof(p,1) .* betaCLM4_5;
% %gs0 is called gw in SoilFLuxesMultiLayer !!!

%REMOVED SEBAS
% % c02 = vegetation.soilvar.thk(c) ./ (vegetation.soilvar.z(c)-vegetation.soilvar.zi(c)); %snl(c))) ; %! Soil heat flux term (W/m2/K
% c02 = vegetation.soilvar.thk_topsurfacecell(c) ./ (vegetation.soilvar.dz_topsurfacecell); %snl(c))) ; %! Soil heat flux term (W/m2/K)
% % c01 = -c02 * vegetation.soilvar.t_soisno(c); %(c,1+1);     %! Soil heat flux term (W/m2)
% c01 = -c02 * vegetation.soilvar.t_top_surfacecell; %(c,1+1);     %! Soil heat flux term (W/m2)

% Coefficients for ground temperature

% den = vegetation.mlcanopyinst.cpair(p) .* vegetation.mlcanopyinst.ga_prof(p,1) + lambda .* vegetation.mlcanopyinst.rhg(p) .* gs0 .* dqsat0 + c02; %p,0
% alpha0 = vegetation.mlcanopyinst.cpair(p) .* vegetation.mlcanopyinst.ga_prof(p,1) ./ den; %p,0
% beta0 = lambda * gs0 ./ den;
% delta0 = (vegetation.mlcanopyinst.rnsoi(p) - lambda .* vegetation.mlcanopyinst.rhg(p) .* gs0 .* (qsat0 - dqsat0 .* vegetation.mlcanopyinst.tair_old(p,1)) - c01) ./ den; %p,0
%end REMOVE SEBAS


% ---------------------------------------------------------------------
% alpha, beta, delta coefficients for leaf temperature:
%
% Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
% Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
%
% and leaf terms needed for scalar conservations equations. These
% are defined for all layers but the leaf terms only exist for
% the layers with leaves.
% ---------------------------------------------------------------------

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)

   if (vegetation.canopy.dpai(p,ic) > 0.)

      % Calculate terms for sunlit and shaded leaves

      for il = 1:vegetation.mlcanopyinst.nleaf

         % Leaf conductances (here these are per unit leaf area)

         gleaf_sh(ic,il) = 2. * vegetation.mlcanopyinst.gbh(p,ic,il);
%          gleaf = vegetation.mlcanopyinst.gs(p,ic,il) * vegetation.mlcanopyinst.gbv(p,ic,il) / (vegetation.mlcanopyinst.gs(p,ic,il) + vegetation.mlcanopyinst.gbv(p,ic,il));
%          gleaf_et(ic,il) = gleaf * vegetation.mlcanopyinst.fdry(p,ic) + vegetation.mlcanopyinst.gbv(p,ic,il) * vegetation.mlcanopyinst.fwet(p,ic);
         gleaf_et(ic,il) = vegetation.mlcanopyinst.gs(p,ic,il)*vegetation.mlcanopyinst.gbv(p,ic,il)/(vegetation.mlcanopyinst.gs(p,ic,il)+vegetation.mlcanopyinst.gbv(p,ic,il)) * vegetation.mlcanopyinst.fdry(p,ic) + vegetation.mlcanopyinst.gbv(p,ic,il) * vegetation.mlcanopyinst.fwet(p,ic);

% Heat capacity of leaves

         heatcap(ic,il) = vegetation.mlcanopyinst.cpleaf(p,ic);

         % Available energy: net radiation

         avail_energy(ic,il) = vegetation.mlcanopyinst.rnleaf(p,ic,il);

         % Saturation vapor pressure and derivative for leaf temperature at time n: Pa -> mol/mol

         [esat, desat] = Satvap (vegetation.mlcanopyinst.tveg_old(p,ic,il));
         qsat = esat / vegetation.mlcanopyinst.pref(p);
         dqsat(ic,il) = desat / vegetation.mlcanopyinst.pref(p);

         % Term for linearized vapor pressure at leaf temperature:
         % qsat(tveg) = qsat(tveg_old) + dqsat * (tveg - tveg_old)
         % Here qsat_term contains the terms with tveg_old

         qsat_term(ic,il) = qsat - dqsat(ic,il) * vegetation.mlcanopyinst.tveg_old(p,ic,il);

         % alpha, beta, delta coefficients for leaf temperature

         den = heatcap(ic,il) / vegetation.params.dtime + gleaf_sh(ic,il) * vegetation.mlcanopyinst.cpair(p) + gleaf_et(ic,il) * lambda * dqsat(ic,il);
         alpha(ic,il) = gleaf_sh(ic,il) * vegetation.mlcanopyinst.cpair(p) / den;
         beta(ic,il) = gleaf_et(ic,il) * lambda / den;
         delta(ic,il) = avail_energy(ic,il) / den ...
                      - lambda * gleaf_et(ic,il) * qsat_term(ic,il) / den ...
                      + heatcap(ic,il) / vegetation.params.dtime * vegetation.mlcanopyinst.tveg_old(p,ic,il) / den;

         % Now scale flux terms for plant area

         if (il == isun)
            pai = vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
         elseif (il == isha)
            pai = vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
         end

         gleaf_sh(ic,il) = gleaf_sh(ic,il) * pai;
         gleaf_et(ic,il) = gleaf_et(ic,il) * pai;
         heatcap(ic,il) = heatcap(ic,il) * pai;
         avail_energy(ic,il) = avail_energy(ic,il) * pai;

      end

   else

      % Zero out terms

      for il = 1:vegetation.mlcanopyinst.nleaf
         gleaf_sh(ic,il) = 0.;
         gleaf_et(ic,il) = 0.;
         heatcap(ic,il) = 0.;
         avail_energy(ic,il) = 0.;
         qsat(ic,il) = 0.;
         dqsat(ic,il) = 0.;
         qsat_term(ic,il) = 0.;
         alpha(ic,il) = 0.;
         beta(ic,il) = 0.;
         delta(ic,il) = 0.;
      end

   end
end

% ---------------------------------------------------------------------
% a,b,c,d coefficients for air temperature:
% a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
%
% a,b,c,d coefficients for water vapor (mole fraction):
% a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
% ---------------------------------------------------------------------

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)

   % Storage term
                                                                                
   rho_dz_over_dt = vegetation.mlcanopyinst.rhomol(p) * (vegetation.mlcanopyinst.zw(p,ic) - vegetation.mlcanopyinst.zw(p,ic-1)) / vegetation.params.dtime;

   % a1,b11,b12,c1,d1 coefficients for air temperature

   
   a1(ic) = -vegetation.mlcanopyinst.ga_prof(p,ic-1); 
   b11(ic) = rho_dz_over_dt + vegetation.mlcanopyinst.ga_prof(p,ic-1) + vegetation.mlcanopyinst.ga_prof(p,ic) ...
           + gleaf_sh(ic,isun) * (1 - alpha(ic,isun)) ...
           + gleaf_sh(ic,isha) * (1 - alpha(ic,isha));
   b12(ic) = -gleaf_sh(ic,isun) * beta(ic,isun) - gleaf_sh(ic,isha) * beta(ic,isha);
   c1(ic) = -vegetation.mlcanopyinst.ga_prof(p,ic);
   d1(ic) = rho_dz_over_dt * vegetation.mlcanopyinst.tair_old(p,ic) + gleaf_sh(ic,isun) * delta(ic,isun) ...
          + gleaf_sh(ic,isha) * delta(ic,isha);

   % Special case for top layer

   if (ic == vegetation.mlcanopyinst.ncan(p))
      c1(ic) = 0.;
      d1(ic) = d1(ic) + vegetation.mlcanopyinst.ga_prof(p,ic) * vegetation.mlcanopyinst.thref(p);
   end

   % Special case for first canopy layer (i.e., immediately above the ground)

   if (ic == 2)
      a1(ic) = 0.;
%       b11(ic) = b11(ic) - vegetation.mlcanopyinst.ga_prof(p,ic-1) .* alpha0;
%       b12(ic) = b12(ic) - vegetation.mlcanopyinst.ga_prof(p,ic-1) .* beta0;
%       d1(ic) = d1(ic) + vegetation.mlcanopyinst.ga_prof(p,ic-1) .* delta0;
    %SEBAS changed!!!
      b11(ic) = b11(ic) - vegetation.mlcanopyinst.ga_prof(p,ic-1);
      d1(ic) = d1(ic) + vegetation.PARENT_SURFACE.STATVAR.Qh ./  vegetation.mlcanopyinst.cpair(p);
      %d1(ic) = d1(ic) + vegetation.mlcanopyinst.ga_prof(p,ic-1) .* vegetation.mlcanopyinst.tg(p);  %term relating to sensible heat flux
   end
% den = vegetation.mlcanopyinst.cpair(p) .* vegetation.mlcanopyinst.ga_prof(p,1) + lambda .* vegetation.mlcanopyinst.rhg(p) .* gs0 .* dqsat0 + c02; %p,0
% alpha0 = vegetation.mlcanopyinst.cpair(p) .* vegetation.mlcanopyinst.ga_prof(p,1) ./ den; %p,0
% beta0 = lambda * gs0 ./ den;
% delta0 = (vegetation.mlcanopyinst.rnsoi(p) - lambda .* vegetation.mlcanopyinst.rhg(p) .* gs0 .* (qsat0 - dqsat0 .* vegetation.mlcanopyinst.tair_old(p,1)) - c01) ./ den; %p,0

   
   
   % a2,b21,b22,c2,d2 coefficients for water vapor (mole fraction)

%    if (ic == 2)
%       ga_prof_ic_minus_one = gs0;
%    else                                                      
    ga_prof_ic_minus_one = vegetation.mlcanopyinst.ga_prof(p,ic-1);
%    end

   a2(ic) = -ga_prof_ic_minus_one;
   b21(ic) = -gleaf_et(ic,isun) * dqsat(ic,isun) * alpha(ic,isun) ...
             -gleaf_et(ic,isha) * dqsat(ic,isha) * alpha(ic,isha);
   b22(ic) = rho_dz_over_dt + ga_prof_ic_minus_one + vegetation.mlcanopyinst.ga_prof(p,ic) ...
           + gleaf_et(ic,isun) * (1 - dqsat(ic,isun) * beta(ic,isun)) ...
           + gleaf_et(ic,isha) * (1 - dqsat(ic,isha) * beta(ic,isha));
   c2(ic) = -vegetation.mlcanopyinst.ga_prof(p,ic);
   d2(ic) = rho_dz_over_dt * vegetation.mlcanopyinst.eair_old(p,ic) / vegetation.mlcanopyinst.pref(p) ...
          + gleaf_et(ic,isun) * (dqsat(ic,isun) * delta(ic,isun) + qsat_term(ic,isun)) ...
          + gleaf_et(ic,isha) * (dqsat(ic,isha) * delta(ic,isha) + qsat_term(ic,isha));

   % Special case for top layer

   if (ic == vegetation.mlcanopyinst.ncan(p))
      c2(ic) = 0.;
      d2(ic) = d2(ic) + vegetation.mlcanopyinst.ga_prof(p,ic) * vegetation.mlcanopyinst.eref(p)/vegetation.mlcanopyinst.pref(p);
   end


   % Special case for first canopy layer (i.e., immediately above the ground)

   if (ic == 2)
      a2(ic) = 0.;
          %SEBAS changed!!!
%       b21(ic) = b21(ic) - gs0 * vegetation.mlcanopyinst.rhg(p) * dqsat0 * alpha0;
%       b22(ic) = b22(ic) - gs0 * vegetation.mlcanopyinst.rhg(p) * dqsat0 * beta0;
%       d2(ic) = d2(ic) + gs0 * vegetation.mlcanopyinst.rhg(p) * (qsat0 + dqsat0 * (delta0 - vegetation.mlcanopyinst.tair_old(p,1))); %p,0
        b22(ic) = b22(ic) - ga_prof_ic_minus_one;
        d2(ic) = d2(ic) + vegetation.PARENT_SURFACE.STATVAR.Qe./ lambda;
%       d2(ic) = d2(ic) + ga_prof_ic_minus_one * vegetation.mlcanopyinst.eg(p) / vegetation.mlcanopyinst.pref(p);
   end
  
end

% ---------------------------------------------------------------------
% Solve for air temperature and water vapor (mole fraction):
%
% a1(i)*T(i-1) + b11(i)*T(i) + b12(i)*q(i) + c1(i)*T(i+1) = d1(i)
% a2(i)*q(i-1) + b21(i)*T(i) + b22(i)*q(i) + c2(i)*q(i+1) = d2(i)
%
% The solution rewrites these equations so that:
% T(i) = f1(i) - e11(i)*T(i+1) - e12(i)*q(i+1) 
% q(i) = f2(i) - e21(i)*T(i+1) - e22(i)*q(i+1) 
% ---------------------------------------------------------------------

% ic = vegetation.canopy.nbot-1;
e11(1) = 0.;
e12(1) = 0.;
e21(1) = 0.;
e22(1) = 0.;
f1(1) = 0.;
f2(1) = 0.;

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)
    
    % The matrix to invert is:
    %
    % B(i)- A(i)*E(i-1)
    %
    % which is a 2x2 matrix. The
    % elements in the 2x2 matrix are:
    %
    %                     | a b |
    % B(i)- A(i)*E(i-1) = | c d |
    %
    % Calculate these elements (denoted by ainv,binv,
    % cinv,dinv) and the determinant of the matrix.
    
    ainv = b11(ic) - a1(ic) * e11(ic-1);
    binv = b12(ic) - a1(ic) * e12(ic-1);
    cinv = b21(ic) - a2(ic) * e21(ic-1);
    dinv = b22(ic) - a2(ic) * e22(ic-1);
    det = ainv * dinv - binv * cinv;
    
    % E(i) = [B(i) - A(i)*E(i-1)]^(-1) * C(i)
    
    e11(ic) = dinv * c1(ic) / det;
    e12(ic) = -binv * c2(ic) / det;
    e21(ic) = -cinv * c1(ic) / det;
    e22(ic) = ainv * c2(ic) / det;
    
    % F(i) = [B(i) - A(i)*E(i-1)]^(-1) * [D(i) - A(i)*F(i-1)]
    
    f1(ic) =  (dinv*(d1(ic) - a1(ic)*f1(ic-1)) - binv*(d2(ic) - a2(ic)*f2(ic-1))) / det;
    f2(ic) = (-cinv*(d1(ic) - a1(ic)*f1(ic-1)) + ainv*(d2(ic) - a2(ic)*f2(ic-1))) / det;
    
end

% Top layer

ic = vegetation.mlcanopyinst.ncan(p);
vegetation.mlcanopyinst.tair(p,ic) = f1(ic);
vegetation.mlcanopyinst.eair(p,ic) = f2(ic);

% Layers through to bottom of canopy

for ic = vegetation.mlcanopyinst.ncan(p)-1: -1: 2
    vegetation.mlcanopyinst.tair(p,ic) = f1(ic) - e11(ic)*vegetation.mlcanopyinst.tair(p,ic+1) - e12(ic)*vegetation.mlcanopyinst.eair(p,ic+1);
    vegetation.mlcanopyinst.eair(p,ic) = f2(ic) - e21(ic)*vegetation.mlcanopyinst.tair(p,ic+1) - e22(ic)*vegetation.mlcanopyinst.eair(p,ic+1);
end



% ---------------------------------------------------------------------
% Calculate leaf temperature:
% Tlsun(i) = alpha_sun(i)*T(i) + beta_sun(i)*q(i) + delta_sun(i)
% Tlsha(i) = alpha_sha(i)*T(i) + beta_sha(i)*q(i) + delta_sha(i)
% ---------------------------------------------------------------------

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)
   vegetation.mlcanopyinst.tveg(p,ic,isun) = alpha(ic,isun)*vegetation.mlcanopyinst.tair(p,ic) ...
                                   + beta(ic,isun)*vegetation.mlcanopyinst.eair(p,ic) + delta(ic,isun);
   vegetation.mlcanopyinst.tveg(p,ic,isha) = alpha(ic,isha)*vegetation.mlcanopyinst.tair(p,ic) ...
                                   + beta(ic,isha)*vegetation.mlcanopyinst.eair(p,ic) + delta(ic,isha);

   % Special checks for no plant area in layer

   if (vegetation.canopy.dpai(p,ic) > 0.)
      pai = vegetation.flux.fracsun(p,ic) * vegetation.canopy.dpai(p,ic);
      if (pai == 0)
         vegetation.mlcanopyinst.tveg(p,ic,isun) = vegetation.mlcanopyinst.tair(p,ic);
      end
      pai = vegetation.flux.fracsha(p,ic) * vegetation.canopy.dpai(p,ic);
      if (pai == 0)
         vegetation.mlcanopyinst.tveg(p,ic,isha) = vegetation.mlcanopyinst.tair(p,ic);
      end
   else
      vegetation.mlcanopyinst.tveg(p,ic,isun) = vegetation.mlcanopyinst.tair(p,ic);
      vegetation.mlcanopyinst.tveg(p,ic,isha) = vegetation.mlcanopyinst.tair(p,ic);
   end
end


%These should no matter
vegetation.mlcanopyinst.tveg(p,1,isun) = vegetation.mlcanopyinst.tair(p,2); %p,0
vegetation.mlcanopyinst.tveg(p,1,isha) = vegetation.mlcanopyinst.tair(p,2); %p,0
vegetation.mlcanopyinst.tair(p,1) = vegetation.mlcanopyinst.tair(p,2); %p,0
vegetation.mlcanopyinst.eair(p,1) = vegetation.mlcanopyinst.eair(p,2);

% !---------------------------------------------------------------------
% ! Convert water vapor from mol/mol to Pa
% !---------------------------------------------------------------------

    for ic = 1:vegetation.mlcanopyinst.ncan(p)
       vegetation.mlcanopyinst.eair(p,ic) = vegetation.mlcanopyinst.eair(p,ic) * vegetation.mlcanopyinst.pref(p);
    end

%  !---------------------------------------------------------------------
%  ! Calculate ground fluxes
%  !---------------------------------------------------------------------




%  Sensible heat flux
%vegetation.mlcanopyinst.shsoi(p) = vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tg(p) - vegetation.mlcanopyinst.tair(p,2)) * vegetation.mlcanopyinst.ga_prof(p,1); %p,1%p,0
vegetation.mlcanopyinst.shsoi(p) = vegetation.PARENT_SURFACE.STATVAR.Qh;
vegetation.mlcanopyinst.sh0 = vegetation.mlcanopyinst.shsoi(p);
%  Latent heat flux
%vegetation.mlcanopyinst.lhsoi(p) = lambda ./ vegetation.mlcanopyinst.pref(p) .*(vegetation.mlcanopyinst.eg(p) - vegetation.mlcanopyinst.eair(p,2)) * gs0; %p,1
vegetation.mlcanopyinst.lhsoi(p) = vegetation.PARENT_SURFACE.STATVAR.Qe;
vegetation.mlcanopyinst.et0 = vegetation.mlcanopyinst.lhsoi(p) ./ lambda ;
vegetation.mlcanopyinst.etsoi(p) = vegetation.mlcanopyinst.lhsoi(p) ./ lambda;

%not sure what that is
%vegetation.mlcanopyinst.g0 = c01 + c02 * t0;



% %     ! Ground temperature
% 
%     t0 = alpha0 * vegetation.mlcanopyinst.tair(p,2) + beta0 * vegetation.mlcanopyinst.eair(p,2) + delta0; %p,1 %p,1
%     vegetation.mlcanopyinst.tveg(p,1,isun) = t0; %p,0
%     vegetation.mlcanopyinst.tveg(p,1,isha) = t0; %p,0
%     vegetation.mlcanopyinst.tair(p,1) = t0; %p,0
% 
% % ! Water vapor at ground (mol/mol)
% 
%     vegetation.mlcanopyinst.eair(p,1) = vegetation.mlcanopyinst.rhg(p) * (qsat0 + dqsat0 * (t0 - vegetation.mlcanopyinst.tair_old(p,1))); %p,0
%     
% % ! Calculate fluxes
% 
%     vegetation.mlcanopyinst.sh0 = -vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tair(p,2) - t0) * vegetation.mlcanopyinst.ga_prof(p,1); %p,1 %p,0
%     vegetation.mlcanopyinst.et0 = -(vegetation.mlcanopyinst.eair(p,2) - vegetation.mlcanopyinst.eair(p,1)) * gs0; %p,1, %p,0
%     vegetation.mlcanopyinst.g0 = c01 + c02 * t0;

% % !---------------------------------------------------------------------
% % ! Convert water vapor from mol/mol to Pa
% % !---------------------------------------------------------------------
% 
%     for ic = 1:vegetation.mlcanopyinst.ncan(p)
%        vegetation.mlcanopyinst.eair(p,ic) = vegetation.mlcanopyinst.eair(p,ic) * vegetation.mlcanopyinst.pref(p);
%     end
    
% ------------------------------------------------------------
% Vertical sensible heat and water vapor fluxes between layers
% ------------------------------------------------------------

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)-1
   vegetation.mlcanopyinst.shair(p,ic) = -vegetation.mlcanopyinst.cpair(p) .* (vegetation.mlcanopyinst.tair(p,ic+1) - vegetation.mlcanopyinst.tair(p,ic)) .* vegetation.mlcanopyinst.ga_prof(p,ic);
   vegetation.mlcanopyinst.etair(p,ic) = -(vegetation.mlcanopyinst.eair(p,ic+1) - vegetation.mlcanopyinst.eair(p,ic)) ./ vegetation.mlcanopyinst.pref(p) .* vegetation.mlcanopyinst.ga_prof(p,ic);
end

ic = vegetation.mlcanopyinst.ncan(p);
vegetation.mlcanopyinst.shair(p,ic) = -vegetation.mlcanopyinst.cpair(p) .* (vegetation.mlcanopyinst.thref(p) - vegetation.mlcanopyinst.tair(p,ic)) * vegetation.mlcanopyinst.ga_prof(p,ic);
vegetation.mlcanopyinst.etair(p,ic) = -(vegetation.mlcanopyinst.eref(p) - vegetation.mlcanopyinst.eair(p,ic)) ./ vegetation.mlcanopyinst.pref(p) .* vegetation.mlcanopyinst.ga_prof(p,ic);

% --------------------------------------------------------------------------
% Canopy air storage flux (W/m2) and its sensible heat and water vapor terms
% --------------------------------------------------------------------------

for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)
   den = (vegetation.mlcanopyinst.zw(p,ic) - vegetation.mlcanopyinst.zw(p,ic-1)) / vegetation.params.dtime;
   storage_sh(ic) = vegetation.mlcanopyinst.rhomol(p) * vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tair(p,ic) - vegetation.mlcanopyinst.tair_old(p,ic)) *den;
   storage_et(ic) = vegetation.mlcanopyinst.rhomol(p) * (vegetation.mlcanopyinst.eair(p,ic) - vegetation.mlcanopyinst.eair_old(p,ic)) ./ vegetation.mlcanopyinst.pref(p) * den;
   vegetation.mlcanopyinst.stair(p,ic) = storage_sh(ic) + storage_et(ic) * lambda;
end

% !---------------------------------------------------------------------
% ! Canopy conservation checks
% !---------------------------------------------------------------------

% ! Canopy layer fluxes as calculated in this routine. Check energy balance
% 
% for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)
%     shsrc(ic) = 0.;
%     etsrc(ic) = 0.;
%     stveg(ic) = 0.;
%     for il = 1:vegetation.mlcanopyinst.nleaf
%         shsrc(ic) = shsrc(ic) + vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tveg(p,ic,il) - vegetation.mlcanopyinst.tair(p,ic)) * gleaf_sh(ic,il);
%         [esat, desat] = Satvap (vegetation.mlcanopyinst.tveg_old(p,ic,il));
%         etsrc(ic) = etsrc(ic) + (esat + desat * (vegetation.mlcanopyinst.tveg(p,ic,il) - vegetation.mlcanopyinst.tveg_old(p,ic,il)) - vegetation.mlcanopyinst.eair(p,ic)) / vegetation.mlcanopyinst.pref(p) * gleaf_et(ic,il);
%         stveg(ic) = stveg(ic) + heatcap(ic,il) * (vegetation.mlcanopyinst.tveg(p,ic,il) - vegetation.mlcanopyinst.tveg_old(p,ic,il)) / vegetation.params.dtime;
%     end
%     
%     error = avail_energy(ic,isun) + avail_energy(ic,isha) - shsrc(ic) - lambda * etsrc(ic) - stveg(ic);
%     if (abs(error) > 0.001)
%         disp(' ERROR: ScalarProfile: Leaf energy balance error');
%         disp(error);
%     end
% end
% 
% % %     ! Flux conservation at each layer. Note special case for first canopy layer.
% 
% for ic = vegetation.canopy.nbot(p):vegetation.mlcanopyinst.ncan(p)
%     
%     if (ic == 2)
%         error = storage_sh(ic) - (vegetation.mlcanopyinst.sh0 + shsrc(ic) - vegetation.mlcanopyinst.shair(p,ic));
%     else
%         error = storage_sh(ic) - (vegetation.mlcanopyinst.shair(p,ic-1) + shsrc(ic) - vegetation.mlcanopyinst.shair(p,ic));
%     end
%     if (abs(error) > 0.001)
%         disp('ERROR: ScalarProfile: Sensible heat layer conservation error');
%         disp(error);
%     end
%     
%     if (ic == 2)
%         error = storage_et(ic) - (vegetation.mlcanopyinst.et0 + etsrc(ic) - vegetation.mlcanopyinst.etair(p,ic));
%     else
%         error = storage_et(ic) - (vegetation.mlcanopyinst.etair(p,ic-1) + etsrc(ic) - vegetation.mlcanopyinst.etair(p,ic));
%     end
%     error = error * lambda;
%     if (abs(error) > 0.001)
%         disp(' ERROR: ScalarProfile: Latent heat layer conservation error');
%         disp(error);
%     end
%     
% end
% 
% %     ! Flux conservation for canopy sensible heat. This is to check canopy
% %     ! conservation equation (so the sum is to ntop not ncan).
% 
% sum_src = 0.;
% sum_storage = 0.;
% for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
%     sum_src = sum_src + shsrc(ic);
%     sum_storage = sum_storage + storage_sh(ic);
% end
% 
% error = (vegetation.mlcanopyinst.sh0 + sum_src - sum_storage) - vegetation.mlcanopyinst.shair(p,vegetation.canopy.ntop(p));
% if (abs(error) > 0.001)
%     disp(' ERROR: ScalarProfile: Sensible heat canopy conservation error');
%         disp(error);
% end
% %     ! Flux conservation for canopy latent heat. This is to check canopy
% %     ! conservation equation (so the sum is to ntop not ncan).
% 
% sum_src = 0.;
% sum_storage = 0.;
% for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
%     sum_src = sum_src + etsrc(ic);
%     sum_storage = sum_storage + storage_et(ic);
% end
% 
% error = (vegetation.mlcanopyinst.et0 + sum_src - sum_storage) - vegetation.mlcanopyinst.etair(p,vegetation.canopy.ntop(p)); %  ! mol H2O/m2/s
% error = error * lambda;                                  %  ! W/m2
% if (abs(error) > 0.001)
%     disp(' ERROR: ScalarProfile: Latent heat canopy conservation error');
%         disp(error);
% end

% %     !---------------------------------------------------------------------
% %     ! Ground temperature energy balance conservation
% %     !---------------------------------------------------------------------

% error = vegetation.mlcanopyinst.rnsoi(p) - vegetation.mlcanopyinst.sh0 - lambda * vegetation.mlcanopyinst.et0 - vegetation.mlcanopyinst.g0;
% if (abs(error) > 0.001)
%     disp(' ERROR: ScalarProfile: Ground temperature energy balance error');
%         disp(error);
% end

% %     !---------------------------------------------------------------------
% %     ! Calculate soil fluxes
% %     !---------------------------------------------------------------------

%[vegetation] = SoilFluxesMultilayer (vegetation, p);

% call SoilFluxesMultilayer (p, soilstate_inst, temperature_inst, &
%     energyflux_inst, waterflux_inst, mlcanopy_inst)
% err = vegetation.mlcanopyinst.shsoi(p)-vegetation.mlcanopyinst.sh0;
% if (abs(err) > 0.001)
% %     disp(' ERROR: ScalarProfile: Soil sensible heat flux error');
% %         disp(err);
% end
% err = lambda*vegetation.mlcanopyinst.etsoi(p)-lambda*vegetation.mlcanopyinst.et0;
% if (abs(err) > 0.001)
% %     disp(' ERROR: ScalarProfile: Soil latent heat flux error');
% %     disp(err);
% end
% % err = vegetation.mlcanopyinst.gsoi(p)-vegetation.mlcanopyinst.g0;
% % if (abs(err) > 0.001)
% % %     disp(' ERROR: ScalarProfile: Soil heat flux error');
% % %     disp(err);
% % end
% err = vegetation.mlcanopyinst.tg(p)-t0;
% if (abs(err) > 1.e-04) %6)
% %     disp(' ERROR: ScalarProfile: Soil temperature error');
% %     disp(err);
% end

end

% % 

% % 

% % 
% %     !---------------------------------------------------------------------
% %     ! Use LeafTemperatureMod to calculate leaf fluxes for the current
% %     ! air temperature and vapor pressure profiles in the canopy. The
% %     ! flux and leaf temperature calculations there are the same as here, so
% %     ! the answers are the same in both routines. But remember that here
% %     ! the fluxes for sunlit/shaded leaves are multiplied by their leaf area.
% %     !---------------------------------------------------------------------
% % 
% %     if (.not. leaf_temp_iter) then
% % 
% %        do ic = 1, ncan(p)
% %           if (dpai(p,ic) > 0) then
% % 
% %              call LeafTemperature (p, ic, isun, mlcanopy_inst)
% %              call LeafTemperature (p, ic, isha, mlcanopy_inst)
% % 
% %              shsrc_leaf = (shleaf(p,ic,isun) * fracsun(p,ic) + shleaf(p,ic,isha) * fracsha(p,ic)) * dpai(p,ic)
% %              etsrc_leaf = (trleaf(p,ic,isun) + evleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) &
% %                         + (trleaf(p,ic,isha) + evleaf(p,ic,isha)) * fracsha(p,ic) * dpai(p,ic)
% %              stveg_leaf = (stleaf(p,ic,isun) * fracsun(p,ic) + stleaf(p,ic,isha) * fracsha(p,ic)) * dpai(p,ic)
% % 
% %              ! Compare flux calculations between LeafTemperatureMod and here
% % 
% %              if (abs(shsrc(ic)-shsrc_leaf) > 0.001_r8) then
% %                 call endrun (msg=' ERROR: ScalarProfile: Leaf sensible heat flux error')
% %              end if
% % 
% %              if (abs(lambda*etsrc(ic)-lambda*etsrc_leaf) > 0.001_r8) then
% %                 call endrun (msg=' ERROR: ScalarProfile: Leaf latent heat flux error')
% %              end if
% % 
% %              if (abs(stveg(ic)-stveg_leaf) > 0.001_r8) then
% %                 call endrun (msg=' ERROR: ScalarProfile: Leaf heat storage error')
% %              end if
% % 
% %              if (abs(tleaf(p,ic,isun)-tveg(p,ic,isun)) > 1.e-06_r8) then
% %                 call endrun (msg=' ERROR: ScalarProfile: Leaf temperature error (sunlit)')
% %              end if
% % 
% %              if (abs(tleaf(p,ic,isha)-tveg(p,ic,isha)) > 1.e-06_r8) then
% %                 call endrun (msg=' ERROR: ScalarProfile: Leaf temperature error (shaded)')
% %              end if
% % 
% %           end if
% %        end do
% % 
% %     end if
% % 
% % !   go to 100
% % 

% % 
% % 100 continue
% % 
% %     end associate
% %   end subroutine ScalarProfile