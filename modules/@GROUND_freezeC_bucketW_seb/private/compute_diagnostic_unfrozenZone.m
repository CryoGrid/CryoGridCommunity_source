function ground = compute_diagnostic_unfrozenZone(ground)


    %adjust for the unfrozen part of the domain
 
    capacityUnfrozen = ground.CONST.c_w .* ground.STATVAR.water ./ ground.STATVAR.layerThick + ...
        ground.CONST.c_m .* ground.STATVAR.mineral ./ ground.STATVAR.layerThick + ground.CONST.c_o .* ground.STATVAR.organic ./ ground.STATVAR.layerThick;
    conductivityUnfrozen = conductivity2(ground.STATVAR.water ./ ground.STATVAR.layerThick, 0, ground.STATVAR.mineral ./ ground.STATVAR.layerThick, ground.STATVAR.organic ./ ground.STATVAR.layerThick, ground);
        
    ground.STATVAR.heatCapacity = double(ground.STATVAR.T<=0) .* ground.STATVAR.heatCapacity + double(ground.STATVAR.T>0) .* capacityUnfrozen;
    ground.STATVAR.thermCond =  double(ground.STATVAR.T<=0) .* ground.STATVAR.thermCond + double(ground.STATVAR.T>0) .* conductivityUnfrozen;
    
   