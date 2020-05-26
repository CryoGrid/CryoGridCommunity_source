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
    ground.STATVAR.heatCapacity = [];
    ground.STATVAR.matricPotential = [];
    %ground.STATVAR.n =[]; %parameter in Van Genuchten parametrization
    %ground.STATVAR.alpha = [];
    %ground.STATVAR.residualWaterContent = [];
    ground.STATVAR.field_capacity = [];
    ground.STATVAR.naturalPorosity = [];
    ground.STATVAR.soil_type = [];
    ground.STATVAR.excessGroundIce = []; % excess ground ice yes/no
     
    ground.STATVAR.Lstar = [];
    ground.STATVAR.Qh = [];
    ground.STATVAR.Qe = [];
    ground.STATVAR.water_reservoir = [];
    ground.STATVAR.surface_runoff = [];
