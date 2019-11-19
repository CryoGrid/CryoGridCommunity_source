function ground = provide_PARA(ground)


ground.PARA.heatFlux_lb = [];

ground.PARA.dt_max = [] ; %[sec]
ground.PARA.dE_max = []; %[J/m3]

%% NC added
%---- initialize the water body module ------------------------------------
% initializeLAKE
ground.PARA.lake.unfrozenWaterSurface = false;
ground.PARA.lake.residualWater = 0; % the water content stored "mixed" cells of air and water if a water body is present


ground.PARA.fieldCapacity=0.50;           % water holding capacity of the soil - this must be adapted to fit the upperlost layers!!
ground.PARA.evaporationDepth=0.10;        % depth to which evaporation occurs - place on grid cell boundaries
ground.PARA.rootDepth=0.20;               % depth affected by transpiration - place on grid cell boundaries
ground.PARA.ratioET=0.5;                  % 1: only transpiration; 0: only evaporation, values in between must be made dependent on LAI, etc.
ground.PARA.externalWaterFlux=0.0;        % external water flux / drainage in [m/day]
ground.PARA.convectiveDomain=[];          % soil domain where air convection due to buoyancy is possible -> start and end [m] - if empty no convection is possible
ground.PARA.mobileWaterDomain=[0 10.0];   % soil domain where water from excess ice melt is mobile -> start and end [m] - if empty water is not mobile
ground.PARA.relative_maxWater=1.0;        % depth at which a water table will form [m] - above excess water is removed, below it pools up
ground.PARA.hydraulic_conductivity = 1e-5;% subsurface saturated hydraulic conductivity assumed for lateral water fluxes [m/s]
ground.PARA.infiltration_limit_depth=2.0; % maxiumum depth [m] from the surface to which infiltration occurse
ground.PARA.arraySize = 0;

%% energy

ground.PARA.energy.E_soil_sens = 0;
ground.PARA.energy.E_soil_lat = 0;
ground.PARA.energy.E_soil = 0;

ground.PARA.energy.E_snow_sens = 0;
ground.PARA.energy.E_snow_lat = 0;
ground.PARA.energy.E_snow = 0;


% WATER ground
% water content soil domain in [m]
ground.PARA.water.W_soil = 0;%nansum( ground.STATVAR.wc .* ground.STATVAR.conductivity );
% water content snow domain in [m]
ground.PARA.water.W_snow = 0;%nansum( GRID.snow.Snow_i + GRID.snow.Snow_w ) + GRID.snow.SWEinitial;
% accumulated changes per output timestep
% storage
ground.PARA.water.dW_soil = 0;
ground.PARA.water.dW_snow = 0;
% precipitation
ground.PARA.water.dp_rain=0;
ground.PARA.water.dp_snow=0; % SWE
% evapotranspiration and sublimation
ground.PARA.water.de=0;
ground.PARA.water.ds=0;
% runoff
ground.PARA.water.dr_surface=0;
ground.PARA.water.dr_external=0;
ground.PARA.water.dr_snowmelt=0;
ground.PARA.water.dr_excessSnow=0;
ground.PARA.water.dr_lateralSnow=0;
ground.PARA.water.dr_rain=0;  % this is only rain on frozen ground
ground.PARA.water.dr_lateralWater=0;
ground.PARA.water.dr_DarcyReservoir=0; % When worker is connected to a Darcy_reservoir as a boundary condition
ground.PARA.water.dr_lateralExcess=0; % excess water when applying lateral fluxes
% mismatch
ground.PARA.water.dm_lacking=0;


% accumulated changes per output timestep
ground.PARA.energy.dE_soil_sens = 0;
ground.PARA.energy.dE_soil_lat = 0;
ground.PARA.energy.dE_soil = 0;
ground.PARA.energy.dE_snow_sens = 0;
ground.PARA.energy.dE_snow_lat = 0;
ground.PARA.energy.dE_snow = 0;

ground.PARA.surf.albedo  = 0.2;
ground.PARA.surf.epsilon = [];
ground.PARA.surf.z0      = [];
ground.PARA.surf.rs      = 50;

ground.PARA.soil.albedo  = [];
ground.PARA.soil.epsilon = [];
ground.PARA.soil.z0 = [];
ground.PARA.soil.rs = [];

ground.PARA.water.albedo  = [];
ground.PARA.water.epsilon = [];
ground.PARA.water.z0 = [];
ground.PARA.water.rs = [];

ground.PARA.ice.albedo  = [];
ground.PARA.ice.epsilon = [];
ground.PARA.ice.z0 = [];
ground.PARA.ice.rs = [];

ground.PARA.technical.z = 2;
ground.PARA.modules.infiltration = 1;

