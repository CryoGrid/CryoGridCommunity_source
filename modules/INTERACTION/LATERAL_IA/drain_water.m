function top_class = drain_water(top_class,mobile_water,waterflux)
    water_out = 0;
    i=1;
    while water_out < abs(waterflux) && sum(mobile_water) > 0
        i = find(mobile_water > 0,1,'first');
        water_out = water_out + mobile_water(i);
        top_class.STATVAR.water(i) = top_class.STATVAR.water(i) - mobile_water(i);
        top_class.STATVAR.waterIce(i) = top_class.STATVAR.waterIce(i) - mobile_water(i);
        mobile_water(i) = 0;
    end
    
    if water_out > abs(waterflux)
        top_class.STATVAR.water(i) = top_class.STATVAR.water(i) + (water_out - abs(waterflux));
        top_class.STATVAR.waterIce(i) = top_class.STATVAR.waterIce(i) + (water_out - abs(waterflux));
    end

end