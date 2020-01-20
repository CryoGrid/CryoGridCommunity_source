function snow_in = get_snow_deposited(drift_index,snow_drifting,lateral)
    
    factor = drift_index(labindex);
    snow_in.ice         = factor.*snow_drifting.ice./lateral.PARA.area(labindex);
    snow_in.waterIce    = factor.*snow_drifting.waterIce./lateral.PARA.area(labindex);
    snow_in.water       = factor.*snow_drifting.water./lateral.PARA.area(labindex);
    snow_in.layerThick  = factor.*snow_drifting.layerThick./lateral.PARA.area(labindex);
    snow_in.energy      = factor.*snow_drifting.energy./lateral.PARA.area(labindex);
    snow_in.d   = snow_drifting.d;
    snow_in.s   = snow_drifting.s;
    snow_in.gs  = snow_drifting.gs;
    snow_in.time_snowfall = snow_drifting.time_snowfall;
    
end