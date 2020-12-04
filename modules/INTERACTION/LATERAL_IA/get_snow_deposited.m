function snow_in = get_snow_deposited(drift_index,snow_drifting)
    
    factor = drift_index(labindex);
    snow_in.ice         = factor.*snow_drifting.ice;
    snow_in.waterIce    = factor.*snow_drifting.waterIce;
    snow_in.water       = factor.*snow_drifting.water;
    snow_in.layerThick  = factor.*snow_drifting.layerThick;
    snow_in.energy      = factor.*snow_drifting.energy;
    snow_in.d   = snow_drifting.d;
    snow_in.s   = snow_drifting.s;
    snow_in.gs  = snow_drifting.gs;
    snow_in.time_snowfall = snow_drifting.time_snowfall;
    
end