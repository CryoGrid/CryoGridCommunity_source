function top_class = add_water_SNOW(top_class,waterflux) % Not compatible for GROUND yet!
water_in = 0;

porespace = top_class.STATVAR.layerThick - (top_class.STATVAR.ice + top_class.STATVAR.water);
% porespace = porespace .* double(top_class.STATVAR.T >=0); % Infiltration only in unfrozen cells

while water_in < waterflux && sum(porespace) > 0
    i = find(porespace > 0,1,'last');
    top_class.STATVAR.water(i) = top_class.STATVAR.water(i) + porespace(i);
    top_class.STATVAR.waterIce(i) = top_class.STATVAR.waterIce(i) + porespace(i);
    water_in = water_in + porespace(i);
    porespace(i) = 0;
    
end

if water_in > waterflux
    top_class.STATVAR.water(i) = top_class.STATVAR.water(i) - (water_in - waterflux);
    top_class.STATVAR.waterIce(i) = top_class.STATVAR.waterIce(i) - (water_in - waterflux);
end


end