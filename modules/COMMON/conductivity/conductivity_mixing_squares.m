function ground = conductivity_mixing_squares(ground)

water = ground.STATVAR.water./ground.STATVAR.layerThick;
ice = ground.STATVAR.ice ./ ground.STATVAR.layerThick;
mineral = ground.STATVAR.mineral./ground.STATVAR.layerThick;
organic = ground.STATVAR.organic./ground.STATVAR.layerThick;
air = ground.STATVAR.air./ground.STATVAR.layerThick;

ground.STATVAR.thermCond = (water.* ground.CONST.k_w.^0.5 + ice.* ground.CONST.k_i.^0.5 ...
    + mineral.* ground.CONST.k_m.^0.5 + organic.* ground.CONST.k_o.^0.5 + air.* ground.CONST.k_a.^0.5).^2;