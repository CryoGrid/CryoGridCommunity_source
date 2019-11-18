function [vegetation] = SoilFluxesMultilayer (vegetation, p)
%p, soilstate_inst, temperature_inst, & energyflux_inst, waterflux_inst, mlcanopy_inst)
%     !
%     ! !DESCRIPTION:
%     ! Soil surface temperature and energy balance
%     !
%     ! !LOCAL VARIABLES:
%     integer  :: c                             ! Column index for CLM g/l/c/p hierarchy
%     real(r8) :: lambda                        ! Latent heat of vaporization (J/mol)
%     real(r8) :: gws                           ! Soil conductance for water vapor (mol H2O/m2/s)
%     real(r8) :: gw                            ! Total conductance for water vapor (mol H2O/m2/s)
%     real(r8) :: esat                          ! Saturation vapor pressure (Pa)
%     real(r8) :: desat                         ! Derivative of saturation vapor pressure (Pa/K)
%     real(r8) :: qsat                          ! Saturation vapor pressure of air (mol/mol)
%     real(r8) :: dqsat                         ! Temperature derivative of saturation vapor pressure (mol/mol/K)
%     real(r8) :: num1, num2, num3, num4, den   ! Intermediate terms
%     real(r8) :: dshsoi                        ! Temperature derivative of sensible heat flux (W/m2/K)
%     real(r8) :: dlhsoi                        ! Temperature derivative of latent heat flux (W/m2/K)
%     real(r8) :: detsoi                        ! Temperature derivative of evaporation flux (mol H2O/m2/s/K)
%     real(r8) :: err                           ! Surface energy imbalance (W/m2)
%     !---------------------------------------------------------------------
%                                                                 ! *** Input ***
%     swsoi          => mlcanopy_inst%swsoi                  , &  ! Absorbed solar radiation, ground (W/m2)
%     irsoi          => mlcanopy_inst%irsoi                  , &  ! Absorbed longwave radiation, ground (W/m2)
%     tref           => mlcanopy_inst%tref                   , &  ! Air temperature at reference height (K)
%     pref           => mlcanopy_inst%pref                   , &  ! Air pressure at reference height (Pa)
%     rhomol         => mlcanopy_inst%rhomol                 , &  ! Molar density at reference height (mol/m3)
%     cpair          => mlcanopy_inst%cpair                  , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
%     tair           => mlcanopy_inst%tair                   , &  ! Air temperature profile (K)
%     eair           => mlcanopy_inst%eair                   , &  ! Vapor pressure profile (Pa)
%     ga_prof        => mlcanopy_inst%ga_prof                , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
%     tair_old       => mlcanopy_inst%tair_old               , &  ! Air temperature profile for previous timestep (K)
%     z              => col%z                                , &  ! Soil layer depth (m)
%     zi             => col%zi                               , &  ! Soil layer depth at layer interface (m)
%     snl            => col%snl                              , &  ! Number of snow layers
%     soilresis      => soilstate_inst%soilresis_col         , &  ! Soil evaporative resistance (s/m)
%     thk            => soilstate_inst%thk_col               , &  ! Soil layer thermal conductivity (W/m/K)
%     smp_l          => soilstate_inst%smp_l_col             , &  ! Soil layer matric potential (mm)
%     t_soisno       => temperature_inst%t_soisno_col        , &  ! Soil temperature (K)
%                                                                 ! *** Output ***
%     rnsoi          => mlcanopy_inst%rnsoi                  , &  ! Net radiation, ground (W/m2)
%     shsoi          => mlcanopy_inst%shsoi                  , &  ! Sensible heat flux, ground (W/m2)
%     lhsoi          => mlcanopy_inst%lhsoi                  , &  ! Latent heat flux, ground (W/m2)
%     gsoi           => mlcanopy_inst%gsoi                   , &  ! Soil heat flux (W/m2)
%     etsoi          => mlcanopy_inst%etsoi                  , &  ! Water vapor flux, ground (mol H2O/m2/s)
%     tg             => mlcanopy_inst%tg                     , &  ! Soil surface temperature (K)
%     eg             => mlcanopy_inst%eg                     , &  ! Soil surface vapor pressure (Pa)
%     rhg            => mlcanopy_inst%rhg                    , &  ! Relative humidity of airspace at soil surface (fraction)
%     eflx_sh_grnd   => energyflux_inst%eflx_sh_grnd_patch   , &  ! CLM: Sensible heat flux from ground (W/m2)
%     eflx_sh_snow   => energyflux_inst%eflx_sh_snow_patch   , &  ! CLM: Sensible heat flux from snow (W/m2)
%     eflx_sh_h2osfc => energyflux_inst%eflx_sh_h2osfc_patch , &  ! CLM: Sensible heat flux from surface water (W/m2)
%     eflx_sh_soil   => energyflux_inst%eflx_sh_soil_patch   , &  ! CLM: Sensible heat flux from soil (W/m2)
%     qflx_evap_soi  => waterflux_inst%qflx_evap_soi_patch   , &  ! CLM: Soil evaporation (kg/m2/s)
%     qflx_ev_snow   => waterflux_inst%qflx_ev_snow_patch    , &  ! CLM: Evaporation flux from snow (kg/m2/s)
%     qflx_ev_h2osfc => waterflux_inst%qflx_ev_h2osfc_patch  , &  ! CLM: Evaporation flux from h2osfc (kg/m2/s)
%     qflx_ev_soil   => waterflux_inst%qflx_ev_soil_patch    , &  ! CLM: Evaporation flux from soil (kg/m2/s)
%     cgrnds         => energyflux_inst%cgrnds_patch         , &  ! CLM: Deriv. of soil sensible heat flux wrt soil temp (W/m2/K)
%     cgrndl         => energyflux_inst%cgrndl_patch         , &  ! CLM: Deriv. of soil evaporation flux wrt soil temp (kg/m2/s/K)
%     cgrnd          => energyflux_inst%cgrnd_patch            &  ! CLM: Deriv. of soil energy flux wrt to soil temp (W/m2/K)

c = p;
ivis = vegetation.params.vis; % Array index for visible waveband
inir = vegetation.params.nir; % Array index for near-infrared waveband

% dz           = 1; %vegetation.soilvar.dz         ;  % Soil layer thickness (m)
% watsat       = 0.4; %vegetation.soilvar.watsat     ;  % Soil layer volumetric water content at saturation (porosity)
% hksat        = 3.; %soilstate_inst.hksat_col      ;  % Soil layer hydraulic conductivity at saturation (mm H2O/s)
% bsw          = 4; %soilstate_inst.bsw_col        ;  % Soil layer Clapp and Hornberger "b" parameter
% rootfr       = 0.5; %soilstate_inst.rootfr_patch   ;  % Fraction of roots in each soil layer
% h2osoi_vol   = 0.1; %waterstate_inst.h2osoi_vol_col;  % Soil layer volumetric water content (m3/m3)
% h2osoi_ice   = 0; %waterstate_inst.h2osoi_ice_col;  % Soil layer ice lens (kg/m2)


%  Current ground temperature

vegetation.mlcanopyinst.tg(p) = vegetation.mlcanopyinst.tair_old(p,1); %p,0

%  Net radiation

vegetation.mlcanopyinst.rnsoi(p) = vegetation.flux.swsoi(p,ivis) + vegetation.flux.swsoi(p,inir) + vegetation.flux.irsoi(p);

%  Latent heat of vaporization

[lambda] = LatVap(vegetation.mlcanopyinst.tref(p), vegetation);

%  Relative humidity in soil airspace
%  Relative humidity of airspace at soil surface (fraction) = exp(Gravitational acceleration (m/s2) * Molecular mass of water (kg/mol) * Soil layer matric potential (mm) / (universal gas constant
%  [J/K/kmole] * Soil temperature (K))
%  fraction = exp(m/s2) * kg/mol * mm / (J/K/kmole * K)
%  vegetation.mlcanopyinst.rhg(p) = exp(vegetation.physcon.grav *  vegetation.physcon.mmh2o *  smp_l(c,1)*1.e-03 / ( vegetation.physcon.rgasc * vegetation.physcon.t_soisno(c,1)));

molar_volume = vegetation.physcon.rgasc * vegetation.soilvar.t_soisno(c,1) / 101.3; % [m3/mol] 
vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential * molar_volume / (vegetation.physcon.rgasc * vegetation.soilvar.t_soisno(c,1)));

% % rgasc back to old value of  8.3145e+03: --> doesn´t matter, molar volume * e+03 / rgasce+03 = same as without e+03
% molar_volume = 8.3145e+03 * vegetation.soilvar.t_soisno(c,1) / 101.3; % [m3/mol] 
% vegetation.mlcanopyinst.rhg(p) = exp(vegetation.soilvar.soil_water_matric_potential * molar_volume / (8.3145e+03 * vegetation.soilvar.t_soisno(c,1)));

%  Soil conductance to water vapour diffusion

gws = 1. / vegetation.mlcanopyinst.soilresis(c);                      % ! s/m -> m/s
gws = gws * vegetation.mlcanopyinst.rhomol(p);                         %   ! m/s -> mol H2O/m2/s
gw = vegetation.mlcanopyinst.ga_prof(p,1) * gws / (vegetation.mlcanopyinst.ga_prof(p,1) + gws);  % ! total conductance

%  Saturation vapor pressure at ground temperature (Pa -> mol/mol)
[esat, desat] = Satvap (vegetation.mlcanopyinst.tg(p));
qsat = esat / vegetation.mlcanopyinst.pref(p) ;
dqsat = desat / vegetation.mlcanopyinst.pref(p);

%  Calculate soil surface temperature
num1 = vegetation.mlcanopyinst.cpair(p) * vegetation.mlcanopyinst.ga_prof(p,1); %p,0
num2 = lambda * gw;
num3 = vegetation.soilvar.thk(c) / (vegetation.soilvar.z(c)-vegetation.soilvar.zi(c)); %no snow layer, so leave out snl(c)+1
num4 = vegetation.mlcanopyinst.rnsoi(p) - num2 * vegetation.mlcanopyinst.rhg(p) * (qsat - dqsat * vegetation.mlcanopyinst.tg(p)) + num3 * vegetation.soilvar.t_soisno(c);
den = num1 + num2 * dqsat * vegetation.mlcanopyinst.rhg(p) + num3;
vegetation.mlcanopyinst.tg(p) = (num1*vegetation.mlcanopyinst.tair(p,2) + num2*vegetation.mlcanopyinst.eair(p,2)/vegetation.mlcanopyinst.pref(p) + num4) / den; %p,1
% vegetation.mlcanopyinst.tg(p)

%  Sensible heat flux

vegetation.mlcanopyinst.shsoi(p) = vegetation.mlcanopyinst.cpair(p) * (vegetation.mlcanopyinst.tg(p) - vegetation.mlcanopyinst.tair(p,2)) * vegetation.mlcanopyinst.ga_prof(p,1); %p,1%p,0

%  Latent heat flux - remember that tair_old(p,0) is tg(p) before the update for time n+1

vegetation.mlcanopyinst.eg(p) = vegetation.mlcanopyinst.rhg(p) * (esat + desat * (vegetation.mlcanopyinst.tg(p) - vegetation.mlcanopyinst.tair_old(p,1))); %p,0
vegetation.mlcanopyinst.lhsoi(p) = lambda / vegetation.mlcanopyinst.pref(p) * (vegetation.mlcanopyinst.eg(p) - vegetation.mlcanopyinst.eair(p,2)) * gw; %p,1

%  Soil heat flux

vegetation.mlcanopyinst.gsoi(p) = vegetation.soilvar.thk(c) * (vegetation.mlcanopyinst.tg(p) - vegetation.soilvar.t_soisno(c)) / (vegetation.soilvar.z(c)-vegetation.soilvar.zi(c));

%  Error check

err = vegetation.mlcanopyinst.rnsoi(p) - vegetation.mlcanopyinst.shsoi(p) - vegetation.mlcanopyinst.lhsoi(p) - vegetation.mlcanopyinst.gsoi(p);
if (abs(err) > 0.001)
    disp(' ERROR: SoilFluxesMultilayerMod: energy balance error');
end

%  Water vapor flux: W/m2 -> mol H2O/m2/s
vegetation.mlcanopyinst.etsoi(p) = vegetation.mlcanopyinst.lhsoi(p) / lambda;

%     ! Needed for CLM soil temperature - Snow fluxes must be set equal to
%     ! fluxes from the soil even if there is no snow
%
%     eflx_sh_grnd(p) = shsoi(p)                ! Sensibe heat flux from ground (W/m2)
%     eflx_sh_snow(p) = shsoi(p)                ! Sensible heat flux from snow (W/m2)
%     eflx_sh_h2osfc(p) = 0._r8                 ! Sensible heat flux from surface water (W/m2)
%     eflx_sh_soil(p) = shsoi(p)                ! Sensible heat flux from soil (W/m2)
%
%     qflx_evap_soi(p) = etsoi(p) * mmh2o       ! Soil evaporation (kg/m2/s)
%     qflx_ev_snow(p) = etsoi(p) * mmh2o        ! Evaporation flux from snow (kg/m2/s)
%     qflx_ev_h2osfc(p) = 0._r8                 ! Evaporation flux from h2osfc (kg/m2/s)
%     qflx_ev_soil(p) = etsoi(p) * mmh2o        ! Evaporation flux from soil (kg/m2/s)
%
%     dshsoi = cpair(p) * ga_prof(p,0)
%     dlhsoi  = lambda / pref(p) * gw * rhg(p) * desat
%     detsoi = dlhsoi / lambda
%
%     cgrnds(p) = dshsoi                        ! Temperature derivative of soil sensible heat flux (W/m2/K)
%     cgrndl(p) = detsoi * mmh2o                ! Temperature deriative of soil evaporation flux (kg/m2/s/K)
%     cgrnd(p) = cgrnds(p) + cgrndl(p)*hvap     ! Temperature derivative of soil energy flux (W/m2/K)

end