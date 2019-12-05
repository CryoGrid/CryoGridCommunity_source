function ground = provide_PARA(ground)
    

    ground.PARA.heatFlux_lb = [];
    
    ground.PARA.dt_max = [] ; %[sec]    
    ground.PARA.dE_max = []; %[J/m3]
    
    ground.PARA.albedo = [];
    ground.PARA.epsilon = [];
    ground.PARA.airT_height = []; %measurement height [m]
    ground.PARA.z0 = []; %roughness length [m]
    ground.PARA.rs = [];
    
    ground.PARA.arraySizeT = []; % number of points in the lookup table
    ground.PARA.kh_mineral = []; %thermal conductivity of mineral
    
    
    