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
