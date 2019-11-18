function ground = finalize_STATVAR(ground)

T = ground.STATVAR.T;
mineral= ground.STATVAR.mineral;
organic= ground.STATVAR.organic;
waterIce= ground.STATVAR.waterIce;
layerThick= ground.STATVAR.layerThick;


% energy = T.*(mineral .* ground.CONST.c_m + organic .* ground.CONST.c_o + double(T>=0).*(waterIce .* ground.CONST.c_w) + ...
%     double(T<0).*(waterIce .* ground.CONST.c_i )) - double(T<0) .* (waterIce) .* ground.CONST.L_f;


ground.STATVAR.waterIce = waterIce .* layerThick; % [m]
ground.STATVAR.mineral = mineral .* layerThick; % [m]
ground.STATVAR.organic = organic .* layerThick; % [m]
%ground.STATVAR.energy = energy .* layerThick;  % [J/m2]

% ground.STATVAR.water = double(T>=0) .* waterIce .* layerThick;  % [m]
% ground.STATVAR.ice = double(T<0) .* waterIce .* layerThick;
ground.STATVAR.air = (1-mineral-organic-waterIce) .* layerThick;  % [m]
% ground = conductivity(ground);

ground.STATVAR.Lstar = -100;
ground.STATVAR.Qh = 0;
ground.STATVAR.Qe = 0;

