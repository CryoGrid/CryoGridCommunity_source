function snow_mixed = get_snow_mixed(lateral)

    snow_mixed.layerThick = 0;
    snow_mixed.ice        = 0;
    snow_mixed.water      = 0;
    snow_mixed.waterIce   = 0;
    snow_mixed.energy     = 0;
    snow_mixed.d          = 0;
    snow_mixed.s          = 0;
    snow_mixed.gs         = 0;
    snow_mixed.time_snowfall = 0;
    
    j = 1;
    while j <= numlabs && lateral.TEMP.ice(j) > 0
        snow_mixed.d          = (snow_mixed.d.*snow_mixed.ice + lateral.TEMP.d(j).*lateral.TEMP.ice(j))./(snow_mixed.ice+lateral.TEMP.ice(j));
        snow_mixed.s          = (snow_mixed.s.*snow_mixed.ice + lateral.TEMP.s(j).*lateral.TEMP.ice(j))./(snow_mixed.ice+lateral.TEMP.ice(j));
        snow_mixed.gs         = (snow_mixed.gs.*snow_mixed.ice + lateral.TEMP.gs(j).*lateral.TEMP.ice(j))./(snow_mixed.ice+lateral.TEMP.ice(j));
        snow_mixed.time_snowfall = (snow_mixed.time_snowfall.*snow_mixed.ice + lateral.TEMP.time_snowfall(j).*lateral.TEMP.ice(j))...
            /(snow_mixed.ice+lateral.TEMP.ice(j));
        snow_mixed.layerThick = snow_mixed.layerThick + lateral.TEMP.layerThick(j);
        snow_mixed.ice        = snow_mixed.ice + lateral.TEMP.ice(j);
        snow_mixed.water      = snow_mixed.water + lateral.TEMP.water(j);
        snow_mixed.waterIce   = snow_mixed.waterIce + lateral.TEMP.waterIce(j);
        snow_mixed.energy     = snow_mixed.energy + lateral.TEMP.energy(j);
        j = j + 1;
    end
end