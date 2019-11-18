function thermCond = conductivity_noSaltDiffusion(ground)

    mineral = ground.STATVAR.mineral;
    porosity = ground.STATVAR.porosity;
    organic = ground.STATVAR.organic;
    a = ground.PARA.a;
    b = ground.PARA.b;
    Tmelt = ground.STATVAR.Tmelt; %in Kelvin
    Tmelt_inDegreeC = Tmelt - 273.15;%convert to °C to compare against temperatures.
    
    T_ub = ground.TEMP.T_ub;

    T = ground.STATVAR.T;

    %conductivity lives on the edges
    layerDepth = [ground.STATVAR.upperPos; ground.STATVAR.upperPos - cumsum(ground.STATVAR.layerThick)];

    %read constants from ground struct

    % W / m K
    k_water = ground.CONST.k_w;
    k_ice = ground.CONST.k_i;
    k_organic = ground.CONST.k_o;
    k_mineral = ground.CONST.k_m;

    %assign thermal conductivity
    %k lives on the edges, i.e. on layerDepth
    %however, we never need the last one
    thermCond = nan(length(ground.STATVAR.layerThick)+1,1);

    %calculate for upper boundary
    i = 1;
    a1 = 1/porosity(i) - a(i)*abs(Tmelt(i))^b(i);

    if (T_ub + T(i))/2 <= Tmelt_inDegreeC(i)
        liqWater = 1/(a1 + a(i)*abs((T_ub + T(i))/2)^b(i));
    else
        liqWater = porosity(i);
    end

    thermCond(i) = (liqWater*sqrt(k_water) + (porosity(i) - liqWater)*sqrt(k_ice) + ...
        mineral(i)*sqrt(k_mineral) + organic(i)*sqrt(k_organic))* ...
        (liqWater*sqrt(k_water) + (porosity(i) - liqWater)*sqrt(k_ice) +  ...
        mineral(i)*sqrt(k_mineral) + organic(i)*sqrt(k_organic));

    %calculate for inner edges
    for i = 2:length(layerDepth)-1

        %take the average of liqWater's of grid cells above and below
        a1=1/porosity(i-1) - a(i-1)*abs(Tmelt(i-1))^b(i-1);
        if T(i-1) <= Tmelt_inDegreeC(i-1)
            liqWater = 1/(a1 + a(i-1)*abs(T(i-1))^b(i-1));
        else
            liqWater = porosity(i-1);
        end

        a1 = 1/porosity(i) - a(i)*abs(Tmelt(i))^b(i);
        if T(i) <= Tmelt(i)
            liqWater = (liqWater + 1/(a1 + a(i)*abs(T(i))^b(i)))/2;
        else
            liqWater = (liqWater+porosity(i))/2;
        end

        thermCond(i) = (liqWater*sqrt(k_water) + ((porosity(i) + porosity(i-1))/2-liqWater)*sqrt(k_ice) + ...
            (mineral(i-1) + mineral(i))/2*sqrt(k_mineral) + ...
            (organic(i-1)+organic(i))/2*sqrt(k_organic))*(liqWater*sqrt(k_water) + ...
            ((porosity(i)+porosity(i-1))/2-liqWater)*sqrt(k_ice) + (mineral(i-1)+mineral(i))/2*sqrt(k_mineral) + ...
            (organic(i-1)+organic(i))/2*sqrt(k_organic));
    end

    %extrapolate for lower edge - this is needed for calculation of
    %interaction
    thermCond(i+1) = thermCond(i);
    
