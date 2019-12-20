function pot_water_fluxes = lateral_Darcy_flux(lateral,pot_water_fluxes,j)
    index = labindex;
    if lateral.PARA.contact_length(j,index) > 0 
        
        if lateral.TEMP.water_table(index) > lateral.TEMP.water_table(j) && lateral.TEMP.mobile_water(index) > 0 % current worker is sending water
            deltaH = lateral.TEMP.water_table(index)-lateral.TEMP.water_table(j);
            contact_height = lateral.TEMP.water_table(index) - max(lateral.TEMP.frost_table(index),lateral.TEMP.frost_table(j));
            distance = sqrt(lateral.PARA.distance(j,index)^2 + (lateral.TEMP.water_table(index)-lateral.TEMP.water_table(j))^2);
            contact_area = lateral.PARA.contact_length(j,index) * contact_height;
            Darcy_flux = lateral.PARA.hydraulic_conductivity(index) * (deltaH/distance) * contact_area; %m^3/s
            pot_water_fluxes(index,j) = Darcy_flux * lateral.PARA.interaction_timestep * 3600 / lateral.PARA.area(j); % m potential water change
            pot_water_fluxes(j,index) = -1* Darcy_flux * lateral.PARA.interaction_timestep * 3600 / lateral.PARA.area(index);
        elseif lateral.TEMP.water_table(index) < lateral.TEMP.water_table(j) && lateral.TEMP.mobile_water(j) > 0 % receiving water
            deltaH = lateral.TEMP.water_table(j)-lateral.TEMP.water_table(index);
            contact_height = lateral.TEMP.water_table(j) - max(lateral.TEMP.frost_table(index),lateral.TEMP.frost_table(j));
            distance = sqrt(lateral.PARA.distance(j,index)^2 + (lateral.TEMP.water_table(index)-lateral.TEMP.water_table(j))^2);
            contact_area = lateral.PARA.contact_length(j,index) * contact_height;
            Darcy_flux = lateral.PARA.hydraulic_conductivity(index) * (deltaH/distance) * contact_area;
            pot_water_fluxes(index,j) = -1* Darcy_flux * lateral.PARA.interaction_timestep * 3600 / lateral.PARA.area(j);
            pot_water_fluxes(j,index) = Darcy_flux * lateral.PARA.interaction_timestep * 3600 / lateral.PARA.area(index);
        else % No water flux
            pot_water_fluxes(index,j) = 0;
            pot_water_fluxes(j,index) = 0;
        end
        
    else % No contact
        pot_water_fluxes(index,j) = 0;
        pot_water_fluxes(j,index) = 0;
    end
    
end