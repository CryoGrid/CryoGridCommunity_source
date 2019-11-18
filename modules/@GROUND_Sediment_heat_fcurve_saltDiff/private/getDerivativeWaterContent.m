function [dliqWater_dT, dliqWater_dsaltConc] = getDerivativeWaterContent(ground)
    %Within this function, T is in Kelvin


    
    T = ground.STATVAR.T + 273.15; %convert from �C to Kelvin
    saltConc = ground.STATVAR.saltConc;
    Tmelt = ground.STATVAR.Tmelt; %in Kelvin
    
    delta = ground.PARA.delta; %timestep for derivative
    R = ground.CONST.R;
    L_f = ground.CONST.L_f;
    Tmelt_free_water = ground.CONST.Tmelt_free_water;

    
    saltConc_unfrozen = -0.5.*(T - Tmelt_free_water).*L_f./R./Tmelt_free_water.^2;

    
    %derivative wrt salt content
    saltConc_left_bound  = min(saltConc - delta, saltConc_unfrozen - delta);
    saltConc_right_bound = min(saltConc + delta, saltConc_unfrozen);
    dliqWater_dsaltConc = double(saltConc < saltConc_unfrozen).* ...
        (watercontent_temperature_salt(ground, {T, saltConc_left_bound}) ...
         - watercontent_temperature_salt(ground, {T, saltConc_right_bound})) ...
        ./(saltConc_left_bound-saltConc_right_bound);
    
    %derivative wrt temperature
    T_left_bound = min(T-delta, Tmelt-delta);
    T_right_bound = min(T+delta, Tmelt);
    dliqWater_dT = double(T < Tmelt).* ...
        (watercontent_temperature_salt(ground, {T_left_bound, saltConc})...
         - watercontent_temperature_salt(ground, {T_right_bound, saltConc})) ...
        ./(T_left_bound-T_right_bound);

end

