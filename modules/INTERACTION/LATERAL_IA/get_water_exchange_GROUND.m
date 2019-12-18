function [exchange, mobile_water] = get_water_exchange_GROUND(top_class)

if top_class.STATVAR.T(1) > 0
    mobile_water = top_class.STATVAR.water - (top_class.STATVAR.layerThick - top_class.STATVAR.ice - top_class.STATVAR.mineral - top_class.STATVAR.organic).*top_class.STATVAR.field_capacity;
    i = find(mobile_water > .001,1,'first'); % min threshld of mobile water for lateral fluxes
    j = i - 1 + find(mobile_water(i:end) <= .001,1,'first'); % last cell with mobile water
    if isempty(j) % Water present throughout the column
        j = length(mobile_water);
    else
        mobile_water(j:end) = 0;
    end
    if i > 1
        mobile_water(1:i) = 0;
    end

    if sum(mobile_water) > 0
        exchange.water_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(i:end));
    else
        exchange.water_table = top_class.STATVAR.lowerPos;
    end
    exchange.frost_table = top_class.STATVAR.lowerPos + sum(top_class.STATVAR.layerThick(j:end));
    exchange.mobile_water = sum(mobile_water(i:j));
else
    exchange.water_table = top_class.STATVAR.lowerPos;
    exchange.frost_table = top_class.STATVAR.lowerPos;
    exchange.mobile_water = 0;
    mobile_water = top_class.STATVAR.water.*0;
end