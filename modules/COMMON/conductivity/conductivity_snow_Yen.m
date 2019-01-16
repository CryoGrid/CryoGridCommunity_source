function snow = conductivity_snow_Yen(snow)


ki = 2.2196 - 0.0062489 .* snow.STATVAR.T + 0.00010154.*snow.STATVAR.T.^2;

snow.STATVAR.thermCond = ki.*(snow.STATVAR.waterIce./snow.STATVAR.layerThick).^1.88;
