function thermCond_snow = thermCond_snow_compiled(ice_snow, layerThick_snow)

snow_density = ice_snow ./ max(1e-20, layerThick_snow) .*920;
thermCond_snow = max(5e-3, 2.3.*(snow_density./1000).^1.88);