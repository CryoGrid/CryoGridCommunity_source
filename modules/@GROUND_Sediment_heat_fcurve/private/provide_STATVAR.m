function ground = provide_STATVAR(ground)

    %additional state variables to the base class
    ground.STATVAR.deltaT =[];
    ground.STATVAR.c_eff = []; %[]
    ground.STATVAR.Tmelt = []; % Kelvin
    ground.STATVAR.soilType = [];
  
    ground.STATVAR.saltConc =[]; %mol
    
    ground.STATVAR.porosity = [];
 