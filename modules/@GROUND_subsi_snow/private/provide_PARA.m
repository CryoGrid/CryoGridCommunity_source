function snow = provide_PARA(snow)
    
    snow.PARA.albedo = [];
    snow.PARA.max_albedo = [];
    snow.PARA.min_albedo = [];
    snow.PARA.epsilon = [];
    snow.PARA.rs = [];
    snow.PARA.z0 = []; %roughness length [m] 
    
    snow.PARA.rho_snow = [];
    snow.PARA.tau_1 = [];
    snow.PARA.tau_a = [];
    snow.PARA.tau_f = [];
    snow.PARA.relative_maxSnow = [];
    snow.PARA.extinction = [];
    snow.PARA.albedo = 0.85; %% NC added for xice module 