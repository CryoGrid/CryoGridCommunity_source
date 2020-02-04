function [snow, snow_out] = get_snow_eroded2(snow,lateral,fraction)
%%%
% snow_out is multiplied by the area 
%%%
    interaction_timestep = lateral.PARA.interaction_timestep.*3600;
    available_snow = snow.TEMP.Si;
    snow_out.layerThick = 0;
    snow_out.ice        = 0;
    snow_out.water      = 0;
    snow_out.waterIce   = 0;
    snow_out.energy     = 0;
    snow_out.d          = 0;
    snow_out.s          = 0;
    snow_out.gs         = 0;
    snow_out.time_snowfall = 0;
    
    if fraction > 1 || fraction < 0
        dege
    end
    
    if sum(double(snow.TEMP.one_over_tau > 0)) > 0
        I = find(snow.TEMP.one_over_tau > 0);
        fractions = min(available_snow(I), snow.TEMP.one_over_tau(I).*interaction_timestep.*fraction.*5);
        weights = snow.STATVAR.ice(I) .* fractions;
        snow_out.layerThick = sum(snow.STATVAR.layerThick(I).* fractions).*lateral.PARA.area(labindex);
        snow_out.ice = sum(snow.STATVAR.ice(I) .* fractions).*lateral.PARA.area(labindex);
        snow_out.water = sum(snow.STATVAR.water(I).*fractions).*lateral.PARA.area(labindex);
        snow_out.waterIce = sum(snow.STATVAR.waterIce(I).*fractions).*lateral.PARA.area(labindex);
        snow_out.energy = sum(snow.STATVAR.energy(I).*fractions).*lateral.PARA.area(labindex);
        snow_out.d = sum(snow.STATVAR.d(I) .* weights)./sum(weights);
        snow_out.s = sum(snow.STATVAR.s(I) .* weights)./sum(weights);
        snow_out.gs = sum(snow.STATVAR.gs(I) .* weights)./sum(weights);
        snow_out.time_snowfall = sum(snow.STATVAR.time_snowfall(I).*weights)./sum(weights);
        
        snow.STATVAR.layerThick(I) = snow.STATVAR.layerThick(I).*(1-fractions);
        snow.STATVAR.ice(I) = snow.STATVAR.ice(I).*(1-fractions);
        snow.STATVAR.water(I) = snow.STATVAR.water(I).*(1-fractions);
        snow.STATVAR.waterIce(I) = snow.STATVAR.waterIce(I).*(1-fractions);
        snow.STATVAR.energy(I) = snow.STATVAR.energy(I).*(1-fractions);
    end

end