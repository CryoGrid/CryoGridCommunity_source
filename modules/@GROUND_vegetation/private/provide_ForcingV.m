function ground = provide_ForcingV(ground)


    %forcing variables need for snow and ground upper boundary
    ground.ForcingV.TEMP.tair = [];
    ground.ForcingV.TEMP.wind = [];
    ground.ForcingV.TEMP.Qh = []; % [m]

    ground.ForcingV.TEMP.Qe = []; % [m]
    ground.ForcingV.TEMP.Sin = []; % [m]
    ground.ForcingV.TEMP.Lin = []; % [m]
    ground.ForcingV.TEMP.p = [];  % [J/m2]
    
    ground.ForcingV.TEMP.snow_reservoir = 0;
    ground.ForcingV.TEMP.snowfall = [];  % [degree C]
    ground.ForcingV.TEMP.rainfall = [];  % [m]
    ground.ForcingV.TEMP.q = [];
    ground.ForcingV.TEMP.t = [];  % [m]
  