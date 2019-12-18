function [exchange, mobile_water] = get_water_exchange_SNOW(top_class)

mobile_water = top_class.STATVAR.water - (top_class.STATVAR.layerThick - top_class.STATVAR.ice).*top_class.PARA.field_capacity;
i = find(mobile_water > .001,1,'first'); % min threshld of mobile water for lateral fluxes
j = find(mobile_water > .001,1,'last'); % last cell with mobile water
mobile_water(mobile_water <= .001) = 0;
if sum(mobile_water) > 0
    exchange.water_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(i+1:end)) + mobile_water(i)/(top_class.STATVAR.layerThick(i) - top_class.STATVAR.ice(i))*top_class.STATVAR.layerThick(i);
else
    exchange.water_table = top_class.STATVAR.lowerPos;
end
exchange.frost_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(j+1:end));
exchange.mobile_water = sum(mobile_water);

end