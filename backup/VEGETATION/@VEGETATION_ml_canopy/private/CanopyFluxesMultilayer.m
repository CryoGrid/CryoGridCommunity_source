% % function [vegetation, p, ic, il] = CanopyFluxesMultilayer(vegetation) %,counter)
% % 
% % % %DESCRIPTION:
% % % Compute fluxes for sunlit and shaded leaves at each level
% % % and for soil surface.
% % % vegetation.physcon.mmh2o = vegetation.physcon.mmh2o;       % Molecular mass of water (kg/mol)
% % % vegetation.physcon.mmdry = vegetation.physcon.mmdry;       % Molecular mass of dry air (kg/mol)
% % % vegetation.physcon.cpd = vegetation.physcon.cpd;           % Specific heat of dry air at constant pressure (J/kg/K)
% % % vegetation.physcon.cpw = vegetation.physcon.cpw;           % Specific heat of water vapor at constant pressure (J/kg/K)
% % % vegetation.physcon.rgasc = vegetation.physcon.rgasc;       % Universal gas constant (J/K/mol)
% % 
% % % vegetation.physcon.num_exposedvegp = vegetation.physcon.num_exposedvegp;       % Number of non-snow-covered veg points in CLM patch filter
% % % vegetation.physcon.filter_exposedvegp = vegetation.physcon.filter_exposedvegp;    % CLM patch filter for non-snow-covered vegetation
% % 
% % %     use clm_time_manager, only : get_nstep, get_step_size
% % %     use SolarRadiationMod, only: SolarRadiation
% % %     use SoilTemperatureMod, only : SoilThermProp
% % %     use CanopyWaterMod, only : CanopyInterception, CanopyEvaporation
% % %     use PlantHydraulicsMod, only: SoilResistance, PlantResistance
% % %     use LeafPhotosynthesisMod, only : PhotosynthesisParam
% % %     use CanopyNitrogenProfileMod, only : CanopyNitrogenProfile
% % %     use LongwaveRadiationMod, only : LongwaveRadiation
% % %     use LeafTemperatureMod, only : LeafHeatCapacity
% % %     use LeafBoundaryLayerMod, only : LeafBoundaryLayer
% % %     use LeafFluxesMod, only : LeafFluxes
% % %     use SoilFluxesMultilayerMod, only : SoilFluxesMultilayer
% % %     use CanopyTurbulenceMod, only : CanopyTurbulence, CanopyTurbulenceDummy
% % %
% % % %ARGUMENTS:
% % %     type(bounds_type), intent(in) :: bounds
% % %     integer, intent(in) :: num_exposedvegp           % Number of non-snow-covered veg points in CLM patch filter
% % %     integer, intent(in) :: filter_exposedvegp(:)     % CLM patch filter for non-snow-covered vegetation
% % %     integer, intent(in) :: num_nolakec               % Number of non-lake points in CLM column filter
% % %     integer, intent(in) :: filter_nolakec(:)         % CLM column filter for non-lake points
% % 
% % %     type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
% % %     type(temperature_type) , intent(inout) :: temperature_inst
% % %     type(waterstate_type)  , intent(inout) :: waterstate_inst
% % %     type(waterflux_type)   , intent(inout) :: waterflux_inst
% % %     type(energyflux_type)  , intent(inout) :: energyflux_inst
% % %     type(frictionvel_type) , intent(inout) :: frictionvel_inst
% % %     type(soilstate_type)   , intent(inout) :: soilstate_inst
% % %     type(surfalb_type)     , intent(inout) :: surfalb_inst
% % %     type(mlcanopy_type)    , intent(inout) :: mlcanopy_inst
% % 
% % % %LOCAL VARIABLES:
% % %     f                                    % Filter index
% % %     p                                    % Patch index for CLM g/l/c/p hierarchy
% % %     c                                    % Column index for CLM g/l/c/p hierarchy
% % %     g                                    % Gridcell index for CLM g/l/c/p hierarchy
% % %     ic                                   % Aboveground layer index
% % %     il                                   % Sunlit (1) or shaded (2) leaf index
% % %     nstep                                % Current model time step number
% % %     dtime                                % Model time step (s)
% % %     num_sub_steps                        % Number of sub-time steps
% % %     niter                                % Current sub-time step
% % 
% % % cv = ones(1,-nlevsno+1,nlevgrnd); (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   % CLM: soil heat capacity (J/m2/K)
% % % tk = ones(1,-nlevsno+1,nlevgrnd); (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   % CLM: soil thermal conductivity at layer interface (W/m/K)
% % % tk_h2osfc = ones(1);(bounds%begc:bounds%endc)                 % CLM: thermal conductivity of h2osfc (W/m/K)
% % 
% % %---------------------------------------------------------------------
% % % *** Miscellaneous input ***
% % 
% % % vegetation.pftcon.minlwp      =vegetation.pftcon.minlwp              ; % Minimum leaf water potential (MPa)
% % % vegetation.mlcanopyinst.lai         =vegetation.mlcanopyinst.lai          ; % Leaf area index of canopy (m2/m2)
% % % vegetation.mlcanopyinst.sai         =vegetation.mlcanopyinst.sai          ; % Stem area index of canopy (m2/m2)
% % % vegetation.mlcanopyinst.ncan        =vegetation.mlcanopyinst.ncan         ; % Number of aboveground layers
% % % vegetation.canopy.nbot         =vegetation.canopy.nbot         ; % Index for bottom leaf layer
% % % vegetation.canopy.ntop        =vegetation.canopy.ntop         ; % Index for top leaf layer
% % % vegetation.mlcanopyinst.dpai        =vegetation.mlcanopyinst.dpai         ; % Layer plant area index (m2/m2)
% % 
% % % tacclim     =vegetation.mlcanopyinst.tacclim      ; % Average air temperature for acclimation (K)
% % % btran       =vegetation.mlcanopyinst.btran        ; % Ball-Berry soil wetness factor (-)
% % 
% % % vegetation.flux.swleaf       =vegetation.flux.swleaf       ; % Leaf absorbed solar radiation (W/m2 leaf)
% % 
% % % *** Atmospheric forcing data ***
% % % vegetation.mlcanopyinst.zref        =vegetation.mlcanopyinst.zref         ; % Reference height (m)
% % % vegetation.mlcanopyinst.zref_old    =vegetation.mlcanopyinst.zref_old     ; % Reference height for previous timestep (m)
% % 
% % % vegetation.mlcanopyinst.tref         =vegetation.mlcanopyinst.tref         ; % Air temperature at reference height (K)
% % % vegetation.mlcanopyinst.qref        =vegetation.mlcanopyinst.qref         ; % Specific humidity at reference height (kg/kg)
% % % vegetation.mlcanopyinst.rhref       =vegetation.mlcanopyinst.rhref        ; % Relative humidity at reference height (%)
% % % vegetation.mlcanopyinst.uref         =vegetation.mlcanopyinst.uref         ; % Wind speed at reference height (m/s)
% % % vegetation.mlcanopyinst.pref         =vegetation.mlcanopyinst.pref         ; % Air pressure at reference height (Pa)
% % % vegetation.atmos.swskyb      =vegetation.atmos.swskyb       ; % Atmospheric direct beam solar radiation (W/m2)
% % % vegetation.atmos.swskyd      =vegetation.atmos.swskyd       ; % Atmospheric diffuse solar radiation (W/m2)
% % 
% % % vegetation.mlcanopyinst.irsky       =vegetation.mlcanopyinst.irsky        ; % Atmospheric longwave radiation (W/m2)
% % % vegetation.mlcanopyinst.qflx_rain   =vegetation.mlcanopyinst.qflx_rain    ; % Rainfall (mm H2O/s = kg H2O/m2/s)
% % % vegetation.mlcanopyinst.qflx_snow    =vegetation.mlcanopyinst.qflx_snow    ; % Snowfall (mm H2O/s = kg H2O/m2/s)
% % % vegetation.mlcanopyinst.co2ref       =vegetation.mlcanopyinst.co2ref       ; % Atmospheric CO2 at reference height (umol/mol)
% % 
% % % vegetation.mlcanopyinst.o2ref       =vegetation.mlcanopyinst.o2ref        ; % Atmospheric O2 at reference height (mmol/mol)
% % % *** Derived atmospheric data ***
% % % vegetation.mlcanopyinst.thref         =vegetation.mlcanopyinst.thref        ; % Atmospheric potential temperature (K)
% % % vegetation.mlcanopyinst.thvref      =vegetation.mlcanopyinst.thvref       ; % Atmospheric virtual potential temperature (K)
% % % vegetation.mlcanopyinst.eref        =vegetation.mlcanopyinst.eref         ; % Vapor pressure at reference height (Pa)
% % % vegetation.mlcanopyinst.rhoair      =vegetation.mlcanopyinst.rhoair       ; % Air density at reference height (kg/m3)
% % % vegetation.mlcanopyinst.rhomol      =vegetation.mlcanopyinst.rhomol       ; % Molar density at reference height (mol/m3)
% % % vegetation.mlcanopyinst.mmair       =vegetation.mlcanopyinst.mmair        ; % Molecular mass of air at reference height (kg/mol)
% % % vegetation.mlcanopyinst.cpair        =vegetation.mlcanopyinst.cpair        ; % Specific heat of air at constant pressure, at reference height (J/mol/K)
% % % *** Calculated ***
% % % vegetation.mlcanopyinst.sumpai      =vegetation.mlcanopyinst.sumpai       ; % Cumulative plant area index (m2/m2)
% % % vegetation.mlcanopyinst.tair        =vegetation.mlcanopyinst.tair         ; % Air temperature profile (K)
% % % vegetation.mlcanopyinst.eair        =vegetation.mlcanopyinst.eair         ; % Vapor pressure profile (Pa)
% % % vegetation.mlcanopyinst.cair        =vegetation.mlcanopyinst.cair         ; % Atmospheric CO2 profile (umol/mol)
% % % vegetation.mlcanopyinst.tveg        =vegetation.mlcanopyinst.tveg         ; % Vegetation temperature profile (K)
% % % vegetation.mlcanopyinst.tair_old    =vegetation.mlcanopyinst.tair_old     ; % Air temperature profile for previous timestep (K)
% % % vegetation.mlcanopyinst.eair_old    =vegetation.mlcanopyinst.eair_old     ; % Vapor pressure profile for previous timestep (Pa)
% % % vegetation.mlcanopyinst.cair_old    =vegetation.mlcanopyinst.cair_old     ; % Atmospheric CO2 profile for previous timestep (umol/mol)
% % % vegetation.mlcanopyinst.tveg_old    =vegetation.mlcanopyinst.tveg_old     ; % Vegetation temperature profile for previous timestep (K)
% % % vegetation.mlcanopyinst.tleaf       =vegetation.mlcanopyinst.tleaf        ; % Leaf temperature (K)
% % % vegetation.mlcanopyinst.tleaf_old   =vegetation.mlcanopyinst.tleaf_old    ; % Leaf temperature for previous timestep (K)
% % % vegetation.mlcanopyinst.shair       =vegetation.mlcanopyinst.shair        ; % Canopy air sensible heat flux (W/m2)
% % % vegetation.mlcanopyinst.etair       =vegetation.mlcanopyinst.etair        ; % Canopy air water vapor flux (mol H2O/m2/s)
% % % vegetation.mlcanopyinst.stair        =vegetation.mlcanopyinst.stair        ; % Canopy air storage heat flux (W/m2)
% % % vegetation.mlcanopyinst.irleaf      =vegetation.mlcanopyinst.irleaf       ; % Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
% % % vegetation.mlcanopyinst.rnleaf      =vegetation.mlcanopyinst.rnleaf       ; % Leaf net radiation (W/m2 leaf)
% % % vegetation.mlcanopyinst.stleaf      =vegetation.mlcanopyinst.stleaf       ; % Leaf storage heat flux (W/m2 leaf)
% % % vegetation.mlcanopyinst.shleaf      =vegetation.mlcanopyinst.shleaf       ; % Leaf sensible heat flux (W/m2 leaf)
% % % vegetation.mlcanopyinst.lhleaf      =vegetation.mlcanopyinst.lhleaf       ; % Leaf latent heat flux (W/m2 leaf)
% % % vegetation.mlcanopyinst.trleaf      =vegetation.mlcanopyinst.trleaf       ; % Leaf transpiration flux (mol H2O/m2 leaf/s)
% % % vegetation.mlcanopyinst.evleaf      =vegetation.mlcanopyinst.evleaf       ; % Leaf evaporation flux (mol H2O/m2 leaf/s)
% % % vegetation.mlcanopyinst.ag          =vegetation.mlcanopyinst.ag           ; % Leaf gross photosynthesis (umol CO2/m2 leaf/s)
% % % vegetation.mlcanopyinst.an          =vegetation.mlcanopyinst.an           ; % Leaf net photosynthesis (umol CO2/m2 leaf/s)
% % % vegetation.mlcanopyinst.gs          =vegetation.mlcanopyinst.gs           ; % Leaf stomatal conductance (mol H2O/m2 leaf/s)
% % % vegetation.flux.fracsun     =vegetation.flux.fracsun      ; % Sunlit fraction of canopy layer
% % % vegetation.flux.fracsha     =vegetation.flux.fracsha      ; % Shaded fraction of canopy layer
% % % vegetation.mlcanopyinst.psil        =vegetation.mlcanopyinst.psil         ; % Leaf water potential (MPa)
% % % vegetation.mlcanopyinst.lwp         =vegetation.mlcanopyinst.lwp          ; % Leaf water potential of canopy layer (MPa)
% % % vegetation.mlcanopyinst.fracminlwp  =vegetation.mlcanopyinst.fracminlwp   ; % Fraction of canopy with lwp < minlwp
% % % vegetation.flux.irsoi       =vegetation.flux.irsoi         ; % Absorbed longwave radiation, ground (W/m2)
% % % vegetation.mlcanopyinst.rnsoi       =vegetation.mlcanopyinst.rnsoi        ; % Net radiation, ground (W/m2)
% % % vegetation.mlcanopyinst.shsoi       =vegetation.mlcanopyinst.shsoi        ; % Sensible heat flux, ground (W/m2)
% % % vegetation.mlcanopyinst.lhsoi       =vegetation.mlcanopyinst.lhsoi        ; % Latent heat flux, ground (W/m2)
% % % vegetation.mlcanopyinst.gsoi        =vegetation.mlcanopyinst.gsoi         ; % Soil heat flux (W/m2)
% % % vegetation.mlcanopyinst.etsoi       =vegetation.mlcanopyinst.etsoi        ; % Water vapor flux, ground (mol H2O/m2/s)
% % % vegetation.mlcanopyinst.ircan       =vegetation.mlcanopyinst.ircan        ; % Upward longwave radiation above canopy (W/m2)
% % % vegetation.mlcanopyinst.ustar       =vegetation.mlcanopyinst.ustar        ; % Friction velocity (m/s)
% % % vegetation.mlcanopyinst.ga_prof     =vegetation.mlcanopyinst.ga_prof      ; % Canopy layer aerodynamic conductance for scalars (mol/m2/s)
% % % vegetation.mlcanopyinst.qflx_prec_intr =vegetation.mlcanopyinst.qflx_prec_intr; % Intercepted precipitation (kg H2O/m2/s)
% % 
% % % *** CLM variables ***
% % %     forc_hgt    =atm2lnd_inst%forc_hgt_grc             ; % CLM: Atmospheric reference height (m)
% % %     forc_u      =atm2lnd_inst%forc_u_grc               ; % CLM: Atmospheric wind speed in east direction (m/s)
% % %     forc_v      =atm2lnd_inst%forc_v_grc               ; % CLM: Atmospheric wind speed in north direction (m/s)
% % %     forc_rh     =atm2lnd_inst%forc_rh_grc              ; % CLM: Atmospheric relative humidity (%)
% % %     forc_pco2   =atm2lnd_inst%forc_pco2_grc            ; % CLM: Atmospheric CO2 partial pressure (Pa)
% % %     forc_po2    =atm2lnd_inst%forc_po2_grc             ; % CLM: Atmospheric O2 partial pressure (Pa)
% % %     forc_solad  =atm2lnd_inst%forc_solad_grc           ; % CLM: Atmospheric direct beam radiation (W/m2)
% % %     forc_solai  =atm2lnd_inst%forc_solai_grc           ; % CLM: Atmospheric diffuse radiation (W/m2)
% % %     forc_t      =atm2lnd_inst%forc_t_forwnscaled_col    ; % CLM: Atmospheric temperature (K)
% % %     forc_q      =atm2lnd_inst%forc_q_forwnscaled_col    ; % CLM: Atmospheric specific humidity (kg/kg)
% % %     forc_pbot   =atm2lnd_inst%forc_pbot_forwnscaled_col ; % CLM: Atmospheric pressure (Pa)
% % %     forc_lwrad  =atm2lnd_inst%forc_lwrad_forwnscaled_col; % CLM: Atmospheric longwave radiation (W/m2)
% % %     forc_rain   =atm2lnd_inst%forc_rain_forwnscaled_col ; % CLM: Rainfall rate (mm/s)
% % %     forc_snow   =atm2lnd_inst%forc_snow_forwnscaled_col ; % CLM: Snowfall rate (mm/s)
% % %     btran_patch =energyflux_inst%btran_patch           ; % CLM: Transpiration wetness factor (0 to 1)
% % %     t_a10_patch =temperature_inst%t_a10_patch             &  % CLM: 10-day running mean of the 2-m temperature (K)
% % 
% % 
% % ivis = vegetation.params.vis; % Array index for visible waveband
% % inir = vegetation.params.nir; % Array index for near-infrared waveband
% % isun = vegetation.params.sun; % Array index for sunlit leaf
% % isha = vegetation.params.sha; % Array index for shaded leaf
% % 
% % 
% % % Get current step (counter) and step size (seconds)
% % % vegetation.canopy.dpai = ones(1, vegetation.mlcanopyinst.ncan);
% % % turb_type = vegetation.physcon.turb_type;
% % 
% % % nstep = 1; %get_nstep();
% % vegetation.params.dtime = 1.0*3600.0; %get_step_size();
% % 
% % % Set time step for sub-stepping of flux calculations and number of sub-steps
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % vegetation.params.dtime_sub = 5.0 * 60.0;
% % num_sub_steps = vegetation.params.dtime/vegetation.params.dtime_sub; %int(dtime / dtime_sub);
% % 
% % %---------------------------------------------------------------------
% % % Soil thermal conductivity and heat capacity - do not need the returned
% % % values. Only need thk(c,snl(c)+1), which is the thermal conductivity of
% % % the first snow/soil layer.
% % %---------------------------------------------------------------------
% % 
% % % should only get data from soil module: Cryo Grid
% % %     call SoilThermProp (bounds, num_nolakec, filter_nolakec, tk(begc:endc,:), cv(begc:endc,:), tk_h2osfc(begc:endc), waterstate_inst, temperature_inst, soilstate_inst)
% % % vegetation.physcon.cv = 1.5e6; %2.5e6; % Soil heat capacity from CryoGrid
% % % vegetation.soilvar.thk = 1.2;
% % 
% % %---------------------------------------------------------------------
% % % Solar radiation transfer through the canopy
% % %---------------------------------------------------------------------
% % 
% % %     call SolarRadiation (bounds, vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, surfalb_inst, mlcanopy_inst)
% % [vegetation] = SolarRadiation(vegetation);
% % %---------------------------------------------------------------------
% % % Plant hydraulics
% % %---------------------------------------------------------------------
% % 
% % % should only get data from soil module: Cryo Grid
% % %     call SoilResistance (vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, soilstate_inst, waterstate_inst, mlcanopy_inst)
% % % if counter < 3
% %     [vegetation] = SoilResistance(vegetation);
% % % end
% % 
% % %     call PlantResistance (num_exposedvegp, vegetation.physcon.filter_exposedvegp, mlcanopy_inst)
% % % if counter < 3
% %     [vegetation] = PlantResistance(vegetation);
% % % end %, p, ic
% % %---------------------------------------------------------------------
% % % Canopy profiles of photosynthetic capacity
% % %---------------------------------------------------------------------
% % 
% % %     call CanopyNitrogenProfile (num_exposedvegp, vegetation.physcon.filter_exposedvegp, mlcanopy_inst)
% % % if counter < 3
% %     [vegetation] = CanopyNitrogenProfile(vegetation);
% % % end %, p, ic
% % 
% % %---------------------------------------------------------------------
% % % Use sub-stepping to calculate fluxes over the full time step
% % %---------------------------------------------------------------------
% % 
% % % Initialize fluxes that are summed over sub-time steps
% % for f = 1:vegetation.canopy.num_exposedvegp
% %     p = vegetation.canopy.filter_exposedvegp(f);
% %     vegetation.mlcanopyinst.qflx_prec_intr_sum(p) = 0 ;
% %     vegetation.mlcanopyinst.ircan_sum(p) = 0.;
% %     vegetation.mlcanopyinst.ustar_sum(p) = 0.;
% %     vegetation.mlcanopyinst.irsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.rnsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.shsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.lhsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.etsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.gsoi_sum(p) = 0.;
% %     vegetation.mlcanopyinst.ga_sum(p,1) = 0.; %Fortran (p,0)
% %     for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran 1:ncan
% %         vegetation.mlcanopyinst.irleaf_sum(p,ic) = 0.;
% %         vegetation.mlcanopyinst.shair_sum(p,ic) = 0.;
% %         vegetation.mlcanopyinst.etair_sum(p,ic) = 0.;
% %         vegetation.mlcanopyinst.stair_sum(p,ic) = 0.;
% %         vegetation.mlcanopyinst.ga_sum(p,ic) = 0.;
% %         for il = 1:vegetation.mlcanopyinst.nleaf
% %             vegetation.mlcanopyinst.rnleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.stleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.shleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.lhleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.trleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.evleaf_sum(p,ic,il) = 0.;
% %             vegetation.mlcanopyinst.an_sum(p,ic,il) = 0 ;
% %             vegetation.mlcanopyinst.ag_sum(p,ic,il) = 0 ;
% %             vegetation.mlcanopyinst.gs_sum(p,ic,il) = 0 ;
% %         end
% %     end
% % end
% % 
% % % Sub-stepping loop
% % for niter = 1:num_sub_steps % this line will not be used when coupled with Cryo Grid
% %     
% %     % Save values for previous timestep
% %     for f = 1:vegetation.canopy.num_exposedvegp
% %         p = vegetation.canopy.filter_exposedvegp(f);
% %         
% %         % 1:ncan (number of aboveground layers), but tleaf/_old are only initialized 1:nlevcan (number leave layers)
% %         
% %         % Vegetation input variables
% %         
% %         % ncan (number of aboveground layers)
% %         % nbot (index for bottom leaf layer)
% %         % ntop (index for top leaf layer)
% %         % nveg (number vegetation layers)
% %         
% %         
% %         
% %         % CLM_VARPAR
% %         
% %         % nlevcan (number leave layers)
% %         % nlevgrnd (number ground layers)
% %         % nlevsoi (number of hydrologically active soil layers)
% % 
% %         % nleaf (number of leaves sunlit/shaded)
% %         % numrad (number of solar radiation bands)
% %         
% %         
% %         
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %             vegetation.mlcanopyinst.tleaf_old(p,ic,isun) = vegetation.mlcanopyinst.tleaf(p,ic,isun);
% %             vegetation.mlcanopyinst.tleaf_old(p,ic,isha) = vegetation.mlcanopyinst.tleaf(p,ic,isha);
% %         end
% %         
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.canopy.ntop
% %             vegetation.mlcanopyinst.lwp(p,ic) = -2.0;
% %             vegetation.mlcanopyinst.h2ocan(p,ic) = 0.;
% %         end
% %         
% %         for ic = 1:vegetation.mlcanopyinst.ncan(p) %Fortran 0:vegetation.mlcanopyinst.ncan(p)
% %             vegetation.mlcanopyinst.tair_old(p,ic) = vegetation.mlcanopyinst.tair(p,ic);
% %             vegetation.mlcanopyinst.eair_old(p,ic) = vegetation.mlcanopyinst.eair(p,ic);
% %             vegetation.mlcanopyinst.cair_old(p,ic) = vegetation.mlcanopyinst.cair(p,ic);
% %             vegetation.mlcanopyinst.tveg_old(p,ic,isun) = vegetation.mlcanopyinst.tveg(p,ic,isun);
% %             vegetation.mlcanopyinst.tveg_old(p,ic,isha) = vegetation.mlcanopyinst.tveg(p,ic,isha);
% %         end
% %     end
% %     
% %     % Canopy interception
% %     %        call CanopyInterception (vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, mlcanopy_inst)
% %     [vegetation] = CanopyInterception(vegetation); %, p, ic, il
% %     
% %     % Longwave radiation transfer through canopy
% %     %        call LongwaveRadiation (bounds, vegetation.physcon.num_exposedvegp, vegetation.physcon.filter_exposedvegp, mlcanopy_inst)
% %     [vegetation] = LongwaveRadiation(vegetation); %, p, ic, il);
% %     % Loop through all patches to calculate fluxes
% %     
% %     for f = 1:vegetation.canopy.num_exposedvegp
% %         p = vegetation.canopy.filter_exposedvegp(f);
% %         
% %         % Photosynthesis parameters
% %         %           call PhotosynthesisParam (p, mlcanopy_inst)
% %         [vegetation] = PhotosynthesisParam (vegetation, p);
% % 
% %         % Leaf fluxes for each canopy layer
% %         %         vegetation.mlcanopyinst.irleaf = ones(1,vegetation.mlcanopyinst.ncan); % why would i have to add this again?`???
% %         
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %             
% %             % Leaf heat capacity
% %             %              call LeafHeatCapacity (p, ic, mlcanopy_inst)
% %             [vegetation] = LeafHeatCapacity(vegetation,p,ic); %, ic  %changed if dpai>0 to a for loop..
% %             % Leaf heat capacity, cpleaf = 784.8. In bonan, 2018, 745
% %             % Jm-2K-1 for grass and deciduous broadleaf. For boreal
% %             % evergreen needleleaf it is 2792 Jm-2K-1.
% %             
% %             % Sunlit leaves: net radiation, boundary layer conductances and fluxes
% %             vegetation.mlcanopyinst.rnleaf(p,ic,isun) = vegetation.flux.swleaf(p,ic,isun,ivis) + vegetation.flux.swleaf(p,ic,isun,inir) + vegetation.mlcanopyinst.irleaf(p,ic);
% %             %              call LeafBoundaryLayer (p, ic, isun, mlcanopy_inst)
% %             [vegetation] = LeafBoundaryLayer (vegetation, p, ic, isun);
% %             %              call LeafFluxes (p, ic, isun, mlcanopy_inst)
% %             [vegetation] = LeafFluxes (vegetation,p,ic,isun);
% %             
% %             
% %             % Shaded leaves: net radiation, boundary layer conductances and fluxes
% %             vegetation.mlcanopyinst.rnleaf(p,ic,isha) = vegetation.flux.swleaf(p,ic,isha,ivis) + vegetation.flux.swleaf(p,ic,isha,inir) + vegetation.mlcanopyinst.irleaf(p,ic);
% %             
% %             %              call LeafBoundaryLayer (p, ic, isha, mlcanopy_inst)
% %             [vegetation] = LeafBoundaryLayer (vegetation,p,ic,isha);
% %             %              call LeafFluxes (p, ic, isun, mlcanopy_inst)
% %             [vegetation] = LeafFluxes (vegetation,p,ic,isha);
% %             
% %         end
% %        
% %         % Soil fluxes
% % 
% %         [vegetation] = SoilFluxesMultilayer(vegetation,p);
% %         
% %         %           call  (p, soilstate_inst, temperature_inst, energyflux_inst, waterflux_inst, mlcanopy_inst)
% %         
% %         % Canopy turbulence, aeorodynamic conductances, and scalar profiles using above-
% %         % and within-canopy coupling with a roughness sublayer (RSL) parameterization
% %         
% %         %         switch(turb_type)
% %         %             case (0)
% %         %                 [vegetation] = CanopyTurbulenceDummy (vegetation, p, ic, il);
% %         %             case (1&2&3&4)
% %         [vegetation] = CanopyTurbulence (vegetation,p,niter);  %dt,
% % 
% %         %         end
% %         
% %         %              case default
% %         %                 call endrun (msg=' ERROR: CanopyFluxesMultilayer: turbulence type not valid')
% %         %           end select
% %         
% %         % Update canopy intercepted water for evaporation and dew
% %         
% %         %           call CanopyEvaporation (p, mlcanopy_inst)
% %         [vegetation] = CanopyEvaporation (vegetation, p); %;, ic, il
% % 
% %         % Canopy water and leaf energy fluxes need to be accumulated over all
% %         % sub-time steps because canopy water is updated every time step.
% %         % All other variables are instantaneous for the final time step.
% %         
% %         vegetation.mlcanopyinst.qflx_prec_intr_sum(p) = vegetation.mlcanopyinst.qflx_prec_intr_sum(p) + vegetation.mlcanopyinst.qflx_prec_intr(p);
% %         vegetation.mlcanopyinst.ircan_sum(p) = vegetation.mlcanopyinst.ircan_sum(p) + vegetation.mlcanopyinst.ircan(p);
% %         vegetation.mlcanopyinst.ustar_sum(p) = vegetation.mlcanopyinst.ustar_sum(p) + vegetation.mlcanopyinst.ustar(p);
% %         vegetation.mlcanopyinst.irsoi_sum(p) = vegetation.mlcanopyinst.irsoi_sum(p) + vegetation.flux.irsoi(p);
% %         vegetation.mlcanopyinst.rnsoi_sum(p) = vegetation.mlcanopyinst.rnsoi_sum(p) + vegetation.mlcanopyinst.rnsoi(p);
% %         vegetation.mlcanopyinst.shsoi_sum(p) = vegetation.mlcanopyinst.shsoi_sum(p) + vegetation.mlcanopyinst.shsoi(p);
% %         vegetation.mlcanopyinst.lhsoi_sum(p) = vegetation.mlcanopyinst.lhsoi_sum(p) + vegetation.mlcanopyinst.lhsoi(p);
% %         vegetation.mlcanopyinst.etsoi_sum(p) = vegetation.mlcanopyinst.etsoi_sum(p) + vegetation.mlcanopyinst.etsoi(p);
% %         vegetation.mlcanopyinst.gsoi_sum(p) = vegetation.mlcanopyinst.gsoi_sum(p) + vegetation.mlcanopyinst.gsoi(p);
% %         vegetation.mlcanopyinst.ga_sum(p,1) = vegetation.mlcanopyinst.ga_sum(p,1) + vegetation.mlcanopyinst.ga_prof(p,1); %p,0 %vegetation.mlcanopyinst.ga_sum(p,0) = vegetation.mlcanopyinst.ga_sum(p,0) + vegetation.mlcanopyinst.ga_prof(p,0);
% %         
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %             vegetation.mlcanopyinst.irleaf_sum(p,ic) = vegetation.mlcanopyinst.irleaf_sum(p,ic) + vegetation.mlcanopyinst.irleaf(p,ic);
% %             vegetation.mlcanopyinst.shair_sum(p,ic) = vegetation.mlcanopyinst.shair_sum(p,ic) + vegetation.mlcanopyinst.shair(p,ic);
% %             vegetation.mlcanopyinst.etair_sum(p,ic) = vegetation.mlcanopyinst.etair_sum(p,ic) + vegetation.mlcanopyinst.etair(p,ic);
% %             vegetation.mlcanopyinst.stair_sum(p,ic) = vegetation.mlcanopyinst.stair_sum(p,ic) + vegetation.mlcanopyinst.stair (p,ic);
% %             vegetation.mlcanopyinst.ga_sum(p,ic) = vegetation.mlcanopyinst.ga_sum(p,ic) + vegetation.mlcanopyinst.ga_prof(p,ic);
% %             for il = 1:vegetation.mlcanopyinst.nleaf
% %                 vegetation.mlcanopyinst.rnleaf_sum(p,ic,il) = vegetation.mlcanopyinst.rnleaf_sum(p,ic,il) + vegetation.mlcanopyinst.rnleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.stleaf_sum(p,ic,il) = vegetation.mlcanopyinst.stleaf_sum(p,ic,il) + vegetation.mlcanopyinst.stleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.shleaf_sum(p,ic,il) = vegetation.mlcanopyinst.shleaf_sum(p,ic,il) + vegetation.mlcanopyinst.shleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.lhleaf_sum(p,ic,il) = vegetation.mlcanopyinst.lhleaf_sum(p,ic,il) + vegetation.mlcanopyinst.lhleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.trleaf_sum(p,ic,il) = vegetation.mlcanopyinst.trleaf_sum(p,ic,il) + vegetation.mlcanopyinst.trleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.evleaf_sum(p,ic,il) = vegetation.mlcanopyinst.evleaf_sum(p,ic,il) + vegetation.mlcanopyinst.evleaf(p,ic,il);
% %                 vegetation.mlcanopyinst.an_sum(p,ic,il) = vegetation.mlcanopyinst.an_sum(p,ic,il) + vegetation.mlcanopyinst.an(p,ic,il);
% %                 vegetation.mlcanopyinst.ag_sum(p,ic,il) = vegetation.mlcanopyinst.ag_sum(p,ic,il) + vegetation.mlcanopyinst.ag(p,ic,il);
% %                 vegetation.mlcanopyinst.gs_sum(p,ic,il) = vegetation.mlcanopyinst.gs_sum(p,ic,il) + vegetation.mlcanopyinst.gs(p,ic,il);
% %             end
% %         end
% %     end % End patch loop
% % end    % End sub-stepping loop
% % 
% % %---------------------------------------------------------------------
% % % Sum leaf and soil fluxes
% % %---------------------------------------------------------------------
% % 
% % for f = 1:vegetation.canopy.num_exposedvegp
% %     p = vegetation.canopy.filter_exposedvegp(f);
% %     
% %     vegetation.mlcanopyinst.qflx_prec_intr(p) = vegetation.mlcanopyinst.qflx_prec_intr_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.ircan(p) = vegetation.mlcanopyinst.ircan_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.ustar(p) = vegetation.mlcanopyinst.ustar_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.flux.irsoi(p) = vegetation.mlcanopyinst.irsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.rnsoi(p) = vegetation.mlcanopyinst.rnsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.shsoi(p) = vegetation.mlcanopyinst.shsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.lhsoi(p) = vegetation.mlcanopyinst.lhsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.etsoi(p) = vegetation.mlcanopyinst.etsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.gsoi(p) = vegetation.mlcanopyinst.gsoi_sum(p) / num_sub_steps; %float(num_sub_steps);
% %     vegetation.mlcanopyinst.ga_prof(p,1) = vegetation.mlcanopyinst.ga_sum(p,1) / num_sub_steps; %float(num_sub_steps);   %p,0          vegetation.mlcanopyinst.ga_prof(p,0) = vegetation.mlcanopyinst.ga_sum(p,0) / num_sub_steps; %float(num_sub_steps);
% %     
% %     for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %         vegetation.mlcanopyinst.irleaf(p,ic) = vegetation.mlcanopyinst.irleaf_sum(p,ic) / num_sub_steps; %float(num_sub_steps);
% %         vegetation.mlcanopyinst.shair(p,ic) = vegetation.mlcanopyinst.shair_sum(p,ic) / num_sub_steps; %float(num_sub_steps);
% %         vegetation.mlcanopyinst.etair(p,ic) = vegetation.mlcanopyinst.etair_sum(p,ic) / num_sub_steps; %float(num_sub_steps);
% %         vegetation.mlcanopyinst.stair (p,ic) = vegetation.mlcanopyinst.stair_sum(p,ic) / num_sub_steps; %float(num_sub_steps);
% %         vegetation.mlcanopyinst.ga_prof(p,ic) = vegetation.mlcanopyinst.ga_sum(p,ic) / num_sub_steps; %float(num_sub_steps);
% %         for il = 1:vegetation.mlcanopyinst.nleaf
% %             vegetation.mlcanopyinst.rnleaf(p,ic,il) = vegetation.mlcanopyinst.rnleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.stleaf(p,ic,il) = vegetation.mlcanopyinst.stleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.shleaf(p,ic,il) = vegetation.mlcanopyinst.shleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.lhleaf(p,ic,il) = vegetation.mlcanopyinst.lhleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.trleaf(p,ic,il) = vegetation.mlcanopyinst.trleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.evleaf(p,ic,il) = vegetation.mlcanopyinst.evleaf_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.an(p,ic,il) = vegetation.mlcanopyinst.an_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.ag(p,ic,il) = vegetation.mlcanopyinst.ag_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %             vegetation.mlcanopyinst.gs(p,ic,il) = vegetation.mlcanopyinst.gs_sum(p,ic,il) / num_sub_steps; %float(num_sub_steps);
% %         end
% %     end
% %     disp([vegetation.mlcanopyinst.rnsoi(p) vegetation.mlcanopyinst.gsoi(p)])
% %     %        call CanopyFluxesSum (p, mlcanopy_inst)
% %     [vegetation] = CanopyFluxesSum (vegetation, p);
% %     
% %     %------------------------------------------------------------------
% %     % Need to merge temperature and leaf water potential for sunlit and
% %     % shaded leaves because sun/shade fractions change over time
% %     %------------------------------------------------------------------
% %     
% %     for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %         if (vegetation.canopy.dpai(p,ic) > 0 )
% %             vegetation.mlcanopyinst.tveg(p,ic,isun) = vegetation.mlcanopyinst.tveg(p,ic,isun) * vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.tveg(p,ic,isha) * vegetation.flux.fracsha(p,ic);
% %             vegetation.mlcanopyinst.tveg(p,ic,isha) = vegetation.mlcanopyinst.tveg(p,ic,isun);
% %         end
% %     end
% %     for ic = vegetation.canopy.nbot (p):vegetation.canopy.ntop(p)
% %         vegetation.mlcanopyinst.tleaf(p,ic,isun) = vegetation.mlcanopyinst.tleaf(p,ic,isun) * vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.tleaf(p,ic,isha) * vegetation.flux.fracsha(p,ic);
% %         vegetation.mlcanopyinst.tleaf(p,ic,isha) = vegetation.mlcanopyinst.tleaf(p,ic,isun);
% %     end
% %     
% %     % Leaf water potential for each canopy layer and fraction of the canopy
% %     % that is water stressed
% %     
% %     % %     switch(vegetation.physcon.gstyp)
% %     % %         case (0&1)
% %     %
% %     
% %     
% %     if (vegetation.params.gstyp == 0)
% %         vegetation.mlcanopyinst.fracminlwp(p) = 0 ;
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %             vegetation.mlcanopyinst.lwp(p,ic) = 0;
% %         end
% %     end
% %     
% %     if (vegetation.params.gstyp == 1)
% %         vegetation.mlcanopyinst.fracminlwp(p) = 0 ;
% %         for ic = vegetation.mlcanopyinst.nbot(p):vegetation.mlcanopyinst.ncan(p) %Fortran ic = 1:ncan
% %             vegetation.mlcanopyinst.lwp(p,ic) = 0;
% %         end
% %     end
% %     
% %     % %         case (2)
% %     if (vegetation.params.gstyp == 2)
% %         vegetation.mlcanopyinst.fracminlwp(p) = 0.;
% %         for ic = vegetation.canopy.nbot(p):vegetation.canopy.ntop(p)
% %             vegetation.mlcanopyinst.lwp(p,ic) = vegetation.mlcanopyinst.psil(p,ic,isun) * vegetation.flux.fracsun(p,ic) + vegetation.mlcanopyinst.psil(p,ic,isha) * vegetation.flux.fracsha(p,ic);
% %             if (vegetation.mlcanopyinst.lwp(p,ic) <= vegetation.leaf.minlwp(p))
% %                 vegetation.mlcanopyinst.fracminlwp(p) = vegetation.mlcanopyinst.fracminlwp(p) + vegetation.canopy.dpai(p,ic);
% %             end
% %         end
% %         
% %         if (vegetation.canopy.lai(p) > 0 )
% %             vegetation.mlcanopyinst.fracminlwp(p) = vegetation.mlcanopyinst.fracminlwp(p) / (vegetation.canopy.lai(p) + vegetation.mlcanopyinst.sai(p));
% %         end
% %     end
% %     % % *** Output! ***
% %     % vegetation.mlcanopyinst.sumpai      = sumpai       ; % Cumulative plant area index (m2/m2)
% %     % vegetation.mlcanopyinst.tair        = tair         ; % Air temperature profile (K)
% %     % vegetation.mlcanopyinst.eair        = eair         ; % Vapor pressure profile (Pa)
% %     % vegetation.mlcanopyinst.cair        = cair         ; % Atmospheric CO2 profile (umol/mol)
% %     % vegetation.mlcanopyinst.tveg        = tveg         ; % Vegetation temperature profile (K)
% %     % vegetation.mlcanopyinst.tair_old    = tair_old     ; % Air temperature profile for previous timestep (K)
% %     % vegetation.mlcanopyinst.eair_old    = eair_old     ; % Vapor pressure profile for previous timestep (Pa)
% %     % vegetation.mlcanopyinst.cair_old    = cair_old     ; % Atmospheric CO2 profile for previous timestep (umol/mol)
% %     % vegetation.mlcanopyinst.tveg_old    = tveg_old     ; % Vegetation temperature profile for previous timestep (K)
% %     % vegetation.mlcanopyinst.tleaf       = tleaf        ; % Leaf temperature (K)
% %     % vegetation.mlcanopyinst.tleaf_old   = tleaf_old    ; % Leaf temperature for previous timestep (K)
% %     % vegetation.mlcanopyinst.shair       = shair        ; % Canopy air sensible heat flux (W/m2)
% %     % vegetation.mlcanopyinst.etair       = etair        ; % Canopy air water vapor flux (mol H2O/m2/s)
% %     % vegetation.mlcanopyinst.stair       = stair        ; % Canopy air storage heat flux (W/m2)
% %     % vegetation.mlcanopyinst.irleaf      = irleaf       ; % Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
% %     % vegetation.mlcanopyinst.rnleaf      = rnleaf       ; % Leaf net radiation (W/m2 leaf)
% %     % vegetation.mlcanopyinst.stleaf      = stleaf       ; % Leaf storage heat flux (W/m2 leaf)
% %     % vegetation.mlcanopyinst.shleaf      = shleaf       ; % Leaf sensible heat flux (W/m2 leaf)
% %     % vegetation.mlcanopyinst.lhleaf      = lhleaf       ; % Leaf latent heat flux (W/m2 leaf)
% %     % vegetation.mlcanopyinst.trleaf      = trleaf       ; % Leaf transpiration flux (mol H2O/m2 leaf/s)
% %     % vegetation.mlcanopyinst.evleaf      = evleaf       ; % Leaf evaporation flux (mol H2O/m2 leaf/s)
% %     % vegetation.mlcanopyinst.ag          = ag           ; % Leaf gross photosynthesis (umol CO2/m2 leaf/s)
% %     % vegetation.mlcanopyinst.an          = an           ; % Leaf net photosynthesis (umol CO2/m2 leaf/s)
% %     % vegetation.mlcanopyinst.gs          = gs           ; % Leaf stomatal conductance (mol H2O/m2 leaf/s)
% %     % vegetation.flux.fracsun             = fracsun      ; % Sunlit fraction of canopy layer
% %     % vegetation.flux.fracsha             = fracsha      ; % Shaded fraction of canopy layer
% %     % vegetation.mlcanopyinst.psil        = psil         ; % Leaf water potential (MPa)
% %     % vegetation.mlcanopyinst.lwp         = lwp          ; % Leaf water potential of canopy layer (MPa)
% %     % vegetation.mlcanopyinst.fracminlwp  = fracminlwp   ; % Fraction of canopy with lwp < minlwp
% %     % vegetation.flux.irsoi               = irsoi         ; % Absorbed longwave radiation, ground (W/m2)
% %     % vegetation.mlcanopyinst.rnsoi       = rnsoi        ; % Net radiation, ground (W/m2)
% %     % vegetation.mlcanopyinst.shsoi       = shsoi        ; % Sensible heat flux, ground (W/m2)
% %     % vegetation.mlcanopyinst.lhsoi       = lhsoi        ; % Latent heat flux, ground (W/m2)
% %     % vegetation.mlcanopyinst.gsoi        = gsoi         ; % Soil heat flux (W/m2)
% %     % vegetation.mlcanopyinst.etsoi       = etsoi        ; % Water vapor flux, ground (mol H2O/m2/s)
% %     % vegetation.mlcanopyinst.ircan       = ircan        ; % Upward longwave radiation above canopy (W/m2)
% %     % vegetation.mlcanopyinst.ustar       = ustar        ; % Friction velocity (m/s)
% %     % vegetation.mlcanopyinst.ga_prof     = ga_prof      ; % Canopy layer aerodynamic conductance for scalars (mol/m2/s)
% %     % vegetation.mlcanopyinst.qflx_prec_intr = qflx_prec_intr; % Intercepted precipitation (kg H2O/m2/s)
% %     
% % end
% % %profile off
