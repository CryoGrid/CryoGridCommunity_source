function thermCond_snow = thermCond_snow_compiled(ice_snow, layerThick_snow)

snow_density = min(1, ice_snow ./ max(1e-20, layerThick_snow));
thermCond_snow = max(5e-2, 2.3.*snow_density.^1.88);
