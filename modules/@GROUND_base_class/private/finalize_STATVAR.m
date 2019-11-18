function ground = finalize_STATVAR(ground)

T = ground.STATVAR.T;
mineral= ground.STATVAR.mineral;
organic= ground.STATVAR.organic;
waterIce= ground.STATVAR.waterIce;
layerThick= ground.STATVAR.layerThick;


energy = T.*(mineral .* ground.CONST.c_m + organic .* ground.CONST.c_o + double(T>=0).*(waterIce .* ground.CONST.c_w) + ...
    double(T<0).*(waterIce .* ground.CONST.c_i )) - double(T<0) .* (waterIce) .* ground.CONST.L_f;


ground.STATVAR.waterIce_tot = waterIce .* layerThick; % [m]
ground.STATVAR.mineral_tot = mineral .* layerThick; % [m]
ground.STATVAR.organic_tot = organic .* layerThick; % [m]
ground.STATVAR.energy_tot = energy .* layerThick;  % [J/m2]

ground.STATVAR.water_tot = double(T>=0) .* waterIce .* layerThick;  % [m]
ground.STATVAR.ice_tot = double(T<0) .* waterIce .* layerThick;
ground.STATVAR.air_tot = (1-mineral-organic-waterIce) .* layerThick;  % [m]

ground.STATVAR.waterIce = waterIce .* (0.*layerThick+1); % [-]
ground.STATVAR.mineral = mineral .* (0.*layerThick+1); % [-]
ground.STATVAR.organic = organic .* (0.*layerThick+1); % [-]
ground.STATVAR.energy = energy .* (0.*layerThick+1);  % [-]

ground.STATVAR.water = double(T>=0) .* waterIce .* (0.*layerThick+1);  % [-]
ground.STATVAR.ice = double(T<0) .* waterIce .* (0.*layerThick+1); %[-]
ground.STATVAR.air = (1-mineral-organic-waterIce) .* (0.*layerThick+1);  % [-]

ground = conductivity(ground);

