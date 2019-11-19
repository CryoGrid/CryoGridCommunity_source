function ground = provide_STATVAR(ground)


%status variables
ground.STATVAR.upperPos = [];
ground.STATVAR.lowerPos = [];
ground.STATVAR.layerThick = []; % [m]

ground.STATVAR.waterIce = []; % [m]
ground.STATVAR.mineral = []; % [m]
ground.STATVAR.organic = []; % [m]
ground.STATVAR.energy = [];  % [J/m2]

ground.STATVAR.T = [];  % [degree C]
ground.STATVAR.water = [];  % [m]
ground.STATVAR.ice = [];
ground.STATVAR.air = [];  % [m]
ground.STATVAR.thermCond = [];

%NC added tracers
ground.STATVAR.field_capacity = [];
ground.STATVAR.midpoints = [];
ground.STATVAR.natPor = [];
ground.STATVAR.excessGroundIce = [];
ground.STATVAR.leaching_rate = [];
ground.STATVAR.total_tracer =[];
ground.STATVAR.total_tracer_fixed  = [];
ground.STATVAR.leaching = [];
ground.STATVAR.gridlength = [];
ground.STATVAR.soilType = [];
ground.STATVAR.cT_frozen = [];
ground.STATVAR.cT_thawed = [];

ground.STATVAR.K_frozen = [];
ground.STATVAR.K_thawed = [];
ground.STATVAR.conductivity = [];%k_c_cTgrid
ground.STATVAR.capacity = []; %c_c_cTgrid = [];
ground.STATVAR.liquidWaterContent = []; %lwc_cTgrid
ground.STATVAR.Q_lateral = [];

ground.STATVAR.k_temp =[];
ground.STATVAR.c_temp =[];
ground.STATVAR.lwc_temp =[];
ground.STATVAR.k_eff = [];

ground.STATVAR.wc = [];  % [m]

