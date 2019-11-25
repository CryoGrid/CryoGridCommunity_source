function water_fluxes = available_lateral_water(lateral,pot_water_fluxes)
    water_fluxes = pot_water_fluxes;
    index = labindex;
    loss = find(pot_water_fluxes(:,index)< 0);
    loss_potential = abs(sum(pot_water_fluxes(loss,index)));
    
    available_water = lateral.TEMP.mobile_water(index);
    
    if available_water >= loss_potential % Enough water available
        scaling = 1;
    else
        scaling = available_water / loss_potential;
    end
    
    i = 1;
    while i <= length(loss)
        water_fluxes(loss(i),index) = scaling*pot_water_fluxes(loss(i),index);
        water_fluxes(index,loss(i)) = scaling*pot_water_fluxes(index,loss(i))*lateral.PARA.area(index)/lateral.PARA.area(loss(i));
        i = i+1;
    end
    
end