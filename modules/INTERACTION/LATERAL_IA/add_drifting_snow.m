function snow = add_drifting_snow(snow,snow_in)

    if snow.STATVAR.ice(1) < snow.PARA.swe_per_cell && snow.STATVAR.ice(1) + snow_in.ice >= snow.PARA.swe_per_cell
        underfill = snow.PARA.swe_per_cell - snow.STATVAR.ice(1);
        fraction = (snow_in.ice - underfill)/snow_in.ice;
        
        snow.STATVAR.ice(1)         = snow.STATVAR.ice(1) + snow_in.ice*(1-fraction);
        snow.STATVAR.water(1)       = snow.STATVAR.water(1) + snow_in.water*(1-fraction);
        snow.STATVAR.waterIce(1)    = snow.STATVAR.waterIce(1) + snow_in.waterIce*(1-fraction);
        snow.STATVAR.layerThick(1)  = snow.STATVAR.layerThick(1) + snow_in.layerThick*(1-fraction);
        snow.STATVAR.energy(1)      = snow.STATVAR.energy(1) + snow_in.energy*(1-fraction);
        snow.STATVAR.d(1)           = (snow.STATVAR.d(1)*snow.STATVAR.ice(1) + snow_in.ice*(1-fraction)*snow_in.d)/(snow_in.ice*(1-fraction)+snow.STATVAR.ice(1));
        snow.STATVAR.s(1)           = (snow.STATVAR.s(1)*snow.STATVAR.ice(1) + snow_in.ice*(1-fraction)*snow_in.s)/(snow_in.ice*(1-fraction)+snow.STATVAR.ice(1));
        snow.STATVAR.gs(1)          = (snow.STATVAR.gs(1)*snow.STATVAR.ice(1) + snow_in.ice*(1-fraction)*snow_in.gs)/(snow_in.ice*(1-fraction)+snow.STATVAR.ice(1));
        snow.STATVAR.time_snowfall(1) = (snow.STATVAR.time_snowfall(1)*snow.STATVAR.ice(1) + snow_in.ice*(1-fraction)*snow_in.time_snowfall)...
            /(snow_in.ice*(1-fraction)+snow.STATVAR.ice(1));
        
        snow_in.ice         = snow_in.ice.*fraction;
        snow_in.water       = snow_in.water.*fraction;
        snow_in.waterIce    = snow_in.waterIce.*fraction;
        snow_in.layerThick  = snow_in.layerThick.*fraction;
        snow_in.energy      = snow_in.energy.*fraction;
    end
    
%     while snow_in.ice(1) >= 1.5*snow.PARA.swe_per_cell
%         fraction = (snow_in.ice(1)-snow.PARA.swe_per_cell)./snow_in.ice(1);
%         
%         snow_in.ice             = [snow_in.ice(1)*fraction; snow_in.ice];
%         snow_in.ice(2)          = snow_in.ice(2)*(1-fraction);
%         snow_in.water           = [snow_in.water(1)*fraction; snow_in.water];
%         snow_in.water(2)        = snow_in.water(2)*(1-fraction);
%         snow_in.waterIce        = [snow_in.waterIce(1)*fraction; snow_in.waterIce];
%         snow_in.waterIce(2)     = snow_in.waterIce(2)*(1-fraction);
%         snow_in.layerThick      = [snow_in.layerThick(1)*fraction; snow_in.layerThick];
%         snow_in.layerThick(2)   = snow_in.layerThick(2)*(1-fraction);
%         snow_in.energy          = [snow_in.energy(1)*fraction; snow_in.energy];
%         snow_in.energy(2)       = snow_in.energy(2)*(1-fraction);
%         snow_in.d               = [snow_in.d(1); snow_in.d];
%         snow_in.s               = [snow_in.s(1); snow_in.s];
%         snow_in.gs              = [snow_in.gs(1); snow_in.gs];
%         snow_in.time_snowfall   = [snow_in.time_snowfall(1); snow_in.time_snowfall];
%     end
    if isempty(snow_in.d) || isnan(snow_in.d)
        dff
    end
    
    if snow_in.waterIce > 0
        snow.STATVAR.ice        = [snow_in.ice; snow.STATVAR.ice];
        snow.STATVAR.water      = [snow_in.water; snow.STATVAR.water];
        snow.STATVAR.waterIce   = [snow_in.waterIce; snow.STATVAR.waterIce];
        snow.STATVAR.layerThick = [snow_in.layerThick; snow.STATVAR.layerThick];
        snow.STATVAR.energy     = [snow_in.energy; snow.STATVAR.energy];
        snow.STATVAR.d          = [snow_in.d; snow.STATVAR.d];
        snow.STATVAR.s          = [snow_in.s; snow.STATVAR.s];
        snow.STATVAR.gs         = [snow_in.gs; snow.STATVAR.gs];
        snow.STATVAR.time_snowfall = [snow_in.time_snowfall; snow.STATVAR.time_snowfall];
        snow.STATVAR.target_density = min(1, snow.STATVAR.ice ./ snow.STATVAR.layerThick);
    end
end
