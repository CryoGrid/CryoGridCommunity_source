function [ground] = initialize_mlcanopyinst(ground, forcing)

ground.STATVAR.vegetation.mlcanopyinst.nleaf = 2;                             % Two leaf types
nleaf = ground.STATVAR.vegetation.mlcanopyinst.nleaf;

ncan = ground.STATVAR.vegetation.mlcanopyinst.ncan;

ground.STATVAR.vegetation.mlcanopyinst.numrad = 2;                            % Number of wavebands
numrad = ground.STATVAR.vegetation.mlcanopyinst.numrad;

ground.STATVAR.vegetation.mlcanopyinst.nlevgrnd = 1;                          % Number of soil layers
nlevgrnd = ground.STATVAR.vegetation.mlcanopyinst.nlevgrnd;

ground.STATVAR.vegetation.mlcanopyinst.pftcon = 1;                             % Plant functional type condition (http://www.cgd.ucar.edu/tss/clm/pfts/pft-physiology.htm)
ground.STATVAR.vegetation.mlcanopyinst.p = 1;
ground.STATVAR.vegetation.mlcanopyinst.f = 1;                                  % Filter exposed ground.STATVAR.vegetation


% ground.STATVAR.vegetation.mlcanopyinst.nlevsno = 1;                            %     ??
% ground.STATVAR.vegetation.mlcanopyinst.lwp = zeros(1,ncan);

ground.STATVAR.vegetation.mlcanopyinst.soilresis = 3361.509423807650; % soilvar.resis(p) = 3361.509423807650;       % Soil evaporative resistance (s/m) --> before 600, why??
ground.STATVAR.vegetation.mlcanopyinst.root_biomass = 250; %500.0;            % (Bonan et al. (2014) Geosci. Model Dev., 7, 2193–2222)                 %Fine root biomass (g biomass / m2)                                ground.STATVAR.vegetation input variables
ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_rain = 0;
ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_grnd_snow = 0;

% ground.STATVAR.vegetation.mlcanopyinst.ic = ones(1); %true(1, ncan);                         %Canopy layer index
% ground.STATVAR.vegetation.mlcanopyinst.il = ones(1); %ones(1, ncan);                         %Sunlit or shaded leaf index

% ground.STATVAR.vegetation.mlcanopyinst.nbot = ones(1);                                       %Index for bottom leaf layer                                        ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.ntop = ones(1);                                      %Index for top leaf layer                                           ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.dlai = ones(1,ncan);                            %Layer leaf area index (m2/m2)                                      ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.dsai = ones(1,ncan);                            %Layer stem area index (m2/m2)                                      ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.dpai = ones(1,ncan);                            %Layer plant area index (m2/m2)                                     ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.sumpai = ones(1,ncan);                          %Cumulative plant area index (m2/m2) [for nlevcan layers]           ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.zs = ones(1,ncan+1);                            %Canopy height for scalar concentration and source (m)              ground.STATVAR.vegetation input variables
% ground.STATVAR.vegetation.mlcanopyinst.zw = zeros(1,ncan-1);                              %Canopy heights at layer interfaces (m)                             ground.STATVAR.vegetation input variables

% Atmospheric input variables
% ground.STATVAR.vegetation.mlcanopyinst.zref = zeros(1);                                 %Reference height (m)                                               Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.zref_old = ground.PARA.zref_old; %15;                                  %Reference height for previous timestep (m)                              Atmospheric input variables

%% SI: Vorher auskommentiert
ground.STATVAR.vegetation.mlcanopyinst.tref = forcing.TEMP.Tair+273.15;                % should be in Kelvin?                       %Air temperature at reference height (K)                            Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.uref = forcing.TEMP.wind;                          %Wind speed at reference height (m/s)                               Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.rhref = 53.871; %ML: this variable is never used!                              %https://github.com/gbonan/bonanmodeling/blob/a10cf764013be58c2def1dbe7c7e52a3213e061e/sp_16_01/sp_16_01.m     %Relative humidity at reference height (                            Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.pref = forcing.TEMP.p;                             %Air pressure at reference height (Pa)                              Atmospheric input variables


ground.STATVAR.vegetation.mlcanopyinst.co2ref = 380;                                   %Atmospheric CO2 at reference height (umol/mol)                     Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.o2ref = 209;                                    %https://en.wikipedia.org/wiki/Atmosphere_of_Earth                  Atmospheric O2 at reference height (mmol/mol)                                      Atmospheric input variables

% [ground.STATVAR.vegetation] = initialize_atmos(ground.STATVAR.vegetation);
ground.STATVAR.vegetation.mlcanopyinst.irsky = forcing.DATA.Lin(1);                          %Atmospheric longwave radiation (W/m2)                              Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.qflx_rain = forcing.DATA.rainfall(1) ./ (24*3600);                 %Rainfall (mm H2O/s = kg H2O/m2/s)                                  Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.qflx_snow = forcing.DATA.snowfall(1) ./ (24*3600);                 %Snowfall (mm H2O/s = kg H2O/m2/s)                                  Atmospheric input variables
ground.STATVAR.vegetation.mlcanopyinst.tacclim = mean(forcing.DATA.Tair(1))+273.15;                 %Average air temperature for acclimation (K)                 Atmospheric input variables

%% SI: Vorher auskommentiert
% Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.eref = 1700;                                      %Vapor pressure at reference height (Pa)                                      Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.qref = 100; %forcing.TEMP.q;                                      %Specific humidity at reference height (kg/kg)                                Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.rhoair = 1.225;                                 %Air density at reference height (kg/m3)                                      Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.rhomol = 1e6;                                    %Molar density at reference height (mol/m3)                                   Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.mmair = 0.02897;                                %https://www.reference.com/science/molecular-weight-air-560257840a10b2e6      Molecular mass of air at reference height (kg/mol)                                      Additional derived input variables
ground.STATVAR.vegetation.mlcanopyinst.cpair = 20;                                     %Specific heat of air at constant pressure, at reference height (J/mol/K)     Additional derived input variables

ground.STATVAR.vegetation.mlcanopyinst.wind = zeros(1,ncan)+forcing.DATA.wind(1);                            %Wind speed profile (m/s)                                                     Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.wind_most = zeros(1,ncan);                       %Wind speed profile from MOST (m/s)                                           Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.tair = zeros(1,ncan)+273.15+forcing.DATA.Tair(1); %ones(1,ncan)*273.15+15;                            %Air temperature profile (K)                                                  Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.tair_most = ones(1,ncan);                       %Air temperature profile from MOST (K)                                        Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.eair = zeros(1,ncan)+2500;                            %Vapor pressure profile (Pa)                                                  Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.cair = zeros(1,ncan)+380;                            %Atmospheric CO2 profile (umol/mol)                                           Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.tveg  = zeros(1,ncan,nleaf)+273.15+forcing.DATA.Tair(1);                    %ground.STATVAR.vegetation temperature profile (K)                                           Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.tair_old = ones(1,ncan)*273.15+10; %15;                        %Air temperature profile for previous timestep (K)                            Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.eair_old = ones(1,ncan)+1700;                        %Vapor pressure profile for previous timestep (Pa)                            Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.cair_old = ones(1,ncan)+380;                        %Atmospheric CO2 profile for previous timestep (umol/mol)                     Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.tveg_old = zeros(1,ncan,nleaf)+273.15+forcing.DATA.Tair(1);                 %ground.STATVAR.vegetation temperature profile for previous timestep (K)                     Canopy layer variables [for nlevcan layers]
%ground.STATVAR.vegetation.mlcanopyinst.fracsun = ones(1,ncan);          %flux         %Sunlit fraction of canopy layer                                        Canopy layer variables [for nlevcan layers]
%ground.STATVAR.vegetation.mlcanopyinst.fracsha = ones(1,ncan);          %flux         %Shaded fraction of canopy layer                                        Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.irleaf = zeros(1,ncan);                          %Leaf absorbed longwave radiation for canopy layer(W/m2 leaf)                 Canopy layer variables [for nlevcan layers]
% ground.STATVAR.vegetation.mlcanopyinst.leafwp = zeros(1,ncan);                    %Leaf water potential of canopy layer (MPa)                             Canopy layer variables [for nlevcan layers]

ground.STATVAR.vegetation.mlcanopyinst.lsc = zeros(1,ncan);                             %Leaf-specific conductance of canopy layer (mmol H2O/m2 leaf/s/MPa)           Canopy layer variables [for nlevcan layers]
% ground.STATVAR.vegetation.mlcanopyinst.h2ocan = zeros(1,ncan);                         %Canopy layer intercepted water (kg H2O/m2)                                  Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.fwet = zeros(1,ncan);                           %Fraction of plant area index that is wet                                    Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.fdry = zeros(1,ncan);                           %Fraction of plant area index that is green and dry                          Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.shair = zeros(1,ncan);                           %Canopy air sensible heat flux (W/m2)                                         Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.etair = zeros(1,ncan);                           %Canopy air water vapor flux (mol H2O/m2/s)                                   Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.stair = zeros(1,ncan);                           %Canopy air storage heat flux (W/m2)                                          Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.sw_prof = zeros(1,ncan,numrad);                  %Canopy layer absorbed solar radiation (W/m2)                                 Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.ir_prof = zeros(1,ncan);                         %Canopy layer absorbed longwave radiation (W/m2)                              Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.rn_prof = zeros(1,ncan);                         %Canopy layer net radiation (W/m2)                                            Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.st_prof = zeros(1,ncan);                         %Canopy layer storage heat flux (W/m2)                                        Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.sh_prof = zeros(1,ncan);                         %Canopy layer sensible heat flux (W/m2)                                       Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.lh_prof = zeros(1,ncan);                         %Canopy layer latent heat flux (W/m2)                                         Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.et_prof = zeros(1,ncan);                         %Canopy layer water vapor flux (mol H2O/m2/s)                                 Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.fc_prof = zeros(1,ncan);                         %Canopy layer CO2 flux (umol CO2/m2/s)                                        Canopy layer variables [for nlevcan layers]
ground.STATVAR.vegetation.mlcanopyinst.ga_prof = zeros(1,ncan);                        %Canopy layer aerodynamic conductance for scalars (mol/m2/s)                 Canopy layer variables [for nlevcan layers]

% Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.tleaf = zeros(1,ncan,nleaf)+273.15+forcing.DATA.Tair(1); %+3;                 %Leaf temperature (K)    Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.tleaf_old = zeros(1,ncan,nleaf)+273.15+forcing.DATA.Tair(1); %+3; %0;             %Leaf temperature for previous timestep (K)                                     Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.rnleaf = zeros(1,ncan,nleaf);                %Leaf net radiation (W/m2 leaf)                                                 Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.stleaf = zeros(1,ncan,nleaf);                %Leaf storage heat flux (W/m2 leaf)                                             Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.shleaf = zeros(1,ncan,nleaf);                %Leaf sensible heat flux (W/m2 leaf)                                            Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.lhleaf = zeros(1,ncan,nleaf);                %Leaf latent heat flux (W/m2 leaf)                                              Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
% ground.STATVAR.vegetation.mlcanopyinst.swleaf = ones(1,ncan,nleaf,numrad);       %Leaf absorbed solar radiation (W/m2 leaf) [for numrad wavebands]               Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.trleaf = zeros(1,ncan,nleaf);                %Leaf transpiration flux (mol H2O/m2 leaf/s)                                    Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.evleaf = zeros(1,ncan,nleaf);                %Leaf evaporation flux (mol H2O/m2 leaf/s)                                      Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.psil = zeros(1,ncan,nleaf);                  %Leaf water potential (MPa)                                                     Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]


%ML: initilize some leaf water potential - needs to be at a better
%place!!! moved from SetUpCanopy.m
ground.STATVAR.vegetation.mlcanopyinst.lwp = zeros(1,ncan,nleaf);


ground.STATVAR.vegetation.mlcanopyinst.gbh = zeros(1,ncan,nleaf);                   %Leaf boundary layer conductance, heat (mol/m2 leaf/s)                          Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.gbv = zeros(1,ncan,nleaf);                   %Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)                       Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.gbc = zeros(1,ncan,nleaf);                  %Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)                       Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
% ground.STATVAR.vegetation.mlcanopyinst.apar = ones(1,ncan,nleaf);    %flux       %Leaf absorbed PAR (umol photon/m2 leaf/s)                                      Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.ac = zeros(1,ncan,nleaf);                    %Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)                 Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.aj = zeros(1,ncan,nleaf);                    %Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)                    Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.ap = zeros(1,ncan,nleaf);                    %Leaf product-limited (C3), CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)    Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.ag = zeros(1,ncan,nleaf);                    %Leaf gross photosynthesis (umol CO2/m2 leaf/s)                                 Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.an = zeros(1,ncan,nleaf);                   %Leaf net photosynthesis (umol CO2/m2 leaf/s)                                   Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.rd = zeros(1,ncan,nleaf);                    %Leaf respiration rate (umol CO2/m2 leaf/s)                                     Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.ci = zeros(1,ncan,nleaf);                    %Leaf intercellular CO2 (umol/mol)                                              Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.cs = zeros(1,ncan,nleaf);                   %Leaf surface CO2 (umol/mol)                                                    Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.hs = zeros(1,ncan,nleaf);                   %Leaf fractional humidity at leaf surface (-)                                   Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.vpd = zeros(1,ncan,nleaf);                   %Leaf vapor pressure deficit (Pa)                                               Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.gs = zeros(1,ncan,nleaf);                   %Leaf stomatal conductance (mol H2O/m2 leaf/s)                                  Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.alphapsn = zeros(1,ncan,nleaf);              %Leaf 13C fractionation factor for photosynthesis (-)                           Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.cpleaf = zeros(1,ncan);                      %Leaf heat capacity (J/m2 leaf/K)                                               Leaf variables [for nlevcan layers and nleaf leaves (sunlit or shaded)]
ground.STATVAR.vegetation.mlcanopyinst.cd = 0.25;                                  % Leaf drag coefficient (dimensionless)

% ground.STATVAR.vegetation.mlcanopyinst.rho = ones(1,numrad);
% ground.STATVAR.vegetation.mlcanopyinst.tau = ones(1,numrad);

% Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.swveg = ones(1,numrad);          %flux                 %Absorbed solar radiation, ground.STATVAR.vegetation (W/m2) [for numrad wavebands]      Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.swvegsun = ones(1,numrad);       %flux                 %Absorbed solar radiation, sunlit canopy (W/m2) [for numrad wavebands]   Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.swvegsha = ones(1,numrad);       %flux                 %Absorbed solar radiation, shaded canopy (W/m2) [for numrad wavebands]   Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.irveg = ones(1);                 %flux                 %Absorbed longwave radiation, ground.STATVAR.vegetation (W/m2)                        Canopy variables (fluxes are per m2 ground area)

% ground.STATVAR.vegetation.mlcanopyinst.swveg = ground.STATVAR.vegetation.flux.swveg;          %flux                 %Absorbed solar radiation, ground.STATVAR.vegetation (W/m2) [for numrad wavebands]      Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.swvegsun = ground.STATVAR.vegetation.flux.swvegsun;       %flux                 %Absorbed solar radiation, sunlit canopy (W/m2) [for numrad wavebands]   Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.swvegsha = ground.STATVAR.vegetation.flux.swvegsha;       %flux                 %Absorbed solar radiation, shaded canopy (W/m2) [for numrad wavebands]   Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.irveg = ground.STATVAR.vegetation.flux.irveg;                 %flux                 %Absorbed longwave radiation, ground.STATVAR.vegetation (W/m2)                        Canopy variables (fluxes are per m2 ground area)

ground.STATVAR.vegetation.mlcanopyinst.irvegsun = zeros(1);                                 %Absorbed longwave radiation, sunlit canopy (W/m2)                     Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.irvegsha = zeros(1);                                 %Absorbed longwave radiation, shaded canopy (W/m2)                     Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.shveg = zeros(1);                              %Sensible heat flux, ground.STATVAR.vegetation (W/m2)                                 Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.shvegsun = zeros(1);                                 %Sensible heat flux, sunlit canopy (W/m2)                              Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.shvegsha = zeros(1);                                 %Sensible heat flux, shaded canopy (W/m2)                              Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.lhveg = zeros(1);                                    %Latent heat flux, ground.STATVAR.vegetation (W/m2)                                   Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.lhvegsun = zeros(1);                                 %Latent heat flux, sunlit canopy (W/m2)                                Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.lhvegsha = zeros(1);                                 %Latent heat flux, shaded canopy (W/m2)                                Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.etveg = zeros(1);                              %Water vapor flux, ground.STATVAR.vegetation (mol H2O/m2/s)                           Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.etvegsun = zeros(1);                                 %Water vapor flux, sunlit canopy (mol H2O/m2/s)                        Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.etvegsha = zeros(1);                                 %Water vapor flux, shaded canopy (mol H2O/m2/s)                        Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.gppveg = zeros(1);                                   %Gross primary production (umol CO2/m2/s)                              Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.gppvegsun = zeros(1);                                %Gross primary production, sunlit canopy (umol CO2/m2/s)               Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.gppvegsha = zeros(1);                                %Gross primary production, shaded canopy (umol CO2/m2/s)               Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.albcan = ones(1,numrad);    %flux                        %Albedo above canopy [for numrad wavebands]                            Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.ircan = zeros(1);                                    %Upward longwave radiation above canopy (W/m2)                         Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.rnet = zeros(1);                                     %Net radiation (W/m2)                                                  Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.stflx = ones(1);                                    %Canopy storage heat flux (W/m2)                                       Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.shflx = zeros(1);                                    %Sensible heat flux (W/m2)                                             Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.lhflx = zeros(1);                                    %Latent heat flux (W/m2)                                               Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.etflx = zeros(1);                                    %Water vapor flux (mol H2O/m2/s)                                       Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.fracminlwp = zeros(1);                               %Fraction of canopy with lwp < minlwp                                  Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.ustar = zeros(1);                                    %Friction velocity (m/s)                                               Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.uforc = ones(1);                                    %Wind speed at reference height including stability effect (m/s)       Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.uaf = zeros(1);                                      %Wind speed at canopy top (m/s)                                        Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.taf = zeros(1);                                      %Air temperature at canopy top (K)                                     Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.qaf = zeros(1);                                      %Specific humidity at canopy top (kg/kg)                               Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.eaf = zeros(1);                                      %Vapor pressure at canopy top (Pa)                                     Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.obu = zeros(1);                                      %Obukhov length (m)                                                    Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.obu_gah = zeros(1);                                  %Obukhov length used for gah (m)                                       Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.obuold = ones(1);                                   %Obukhov length from previous iteration                                Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.nmozsgn = ones(1);                                  %Number of times stability changes sign during iteration               Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.z0mg = zeros(1);                                     %Roughness length of ground (m)                                        Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.thref = zeros(1);                                    %Atmospheric potential temperature (K)                                 Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.thvref = zeros(1);                                   %Atmospheric virtual potential temperature (K)                         Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.gah = zeros(1);                                      %Aerodynamic conductance for a scalar above canopy (mol/m2/s)          Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.PrSc = zeros(1);                                     %Prandtl (Schmidt) number at canopy top                                Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.Lc = zeros(1);                                       %Canopy density length scale (m)                                       Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.zdisp = zeros(1);                                    %Displacement height (m)                                               Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.tstar = zeros(1);                                    %Temperature scale (K)                                                 Canopy variables (fluxes are per m2 ground area)
ground.STATVAR.vegetation.mlcanopyinst.qstar = zeros(1);                                    %Water vapor scale (kg/kg)                                             Canopy variables (fluxes are per m2 ground area)
% ground.STATVAR.vegetation.mlcanopyinst.td = zeros(1,ncan);           % http://www.cgd.ucar.edu/tss/clm/pfts/pft-physiology.htm                       %Exponential transmittance of diffuse radiation through a single leaf layer      Canopy variables (fluxes are per m2 ground area)

% Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.rnsoi = zeros(1);                                    %Net radiation, ground (W/m2)                                           Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.shsoi = zeros(1);                                    %Sensible heat flux, ground (W/m2)                                      Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.lhsoi = zeros(1);                                    %Latent heat flux, ground (W/m2)                                        Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.gsoi = zeros(1);                                     %Soil heat flux (W/m2)                                                  Soil energy balance

% ground.STATVAR.vegetation.mlcanopyinst.swsoi = ones(1,numrad);           %flux                  %Absorbed solar radiation, ground (W/m2) [for numrad wavebands]         Soil energy balance
% ground.STATVAR.vegetation.mlcanopyinst.irsoi = ones(1);                  %flux                  %Absorbed longwave radiation, ground (W/m2)                             Soil energy balance

% ground.STATVAR.vegetation.mlcanopyinst.swsoi = ground.STATVAR.vegetation.flux.swsoi;           %flux                  %Absorbed solar radiation, ground (W/m2) [for numrad wavebands]         Soil energy balance
% ground.STATVAR.vegetation.mlcanopyinst.irsoi = ground.STATVAR.vegetation.flux.irsoi;                  %flux                  %Absorbed longwave radiation, ground (W/m2)                             Soil energy balance

ground.STATVAR.vegetation.mlcanopyinst.etsoi = zeros(1);                                    %Water vapor flux, ground (mol H2O/m2/s)                                Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.tg = ground.STATVAR.vegetation.mlcanopyinst.tair(1);               %Soil surface temperature (K)                                           Soil energy balance
ground.STATVAR.vegetation.mlcanopyinst.tg_snow = ground.STATVAR.vegetation.mlcanopyinst.tair(1);  % Soil moisture variables
ground.STATVAR.vegetation.mlcanopyinst.snow = 0;


%SEBAS: implemented again, compare Lines 141:151 in LeafPhotosynthesis.m
 ground.STATVAR.vegetation.mlcanopyinst.btran = 0.8; %10; %1; %0.1; %0.6;     %soil water transpiration factor (0 to 1) %Ball-Berry soil wetness factor (-)                                     Soil moisture variables    %http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesmBbrowser/html_code/clm/CanopyFluxesMod.F90.html#SHR_CONST_MOD_23
%end change SEBAS

 ground.STATVAR.vegetation.mlcanopyinst.psis = zeros(1);                                     %Weighted soil water potential (MPa)                                    Soil moisture variables
ground.STATVAR.vegetation.mlcanopyinst.rsoil = zeros(1);                                    %Soil hydraulic resistance (MPa.s.m2/mmol H2O)                          Soil moisture variables
ground.STATVAR.vegetation.mlcanopyinst.soil_et_loss = zeros(1,nlevgrnd);                    %Fraction of total transpiration from each soil layer (-)               Soil moisture variables

% Water flux variables
ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_intr = zeros(1);                           %Intercepted precipitation (kg H2O/m2/s)                                Water flux variables

ground.STATVAR.vegetation.mlcanopyinst.qflx_prec_intr_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.ircan_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.ustar_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.irsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.rnsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.shsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.lhsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.etsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.gsoi_sum = zeros(1);
ground.STATVAR.vegetation.mlcanopyinst.ga_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan); %In Fortran: 0:nlevcan
ground.STATVAR.vegetation.mlcanopyinst.shair_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.mlcanopyinst.etair_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.mlcanopyinst.stair_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.mlcanopyinst.irleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan);
ground.STATVAR.vegetation.mlcanopyinst.rnleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.stleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.shleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.lhleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.trleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.evleaf_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.an_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.ag_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
ground.STATVAR.vegetation.mlcanopyinst.gs_sum = zeros(1,ground.STATVAR.vegetation.mlcanopyinst.ncan,ground.STATVAR.vegetation.mlcanopyinst.nleaf);
end
