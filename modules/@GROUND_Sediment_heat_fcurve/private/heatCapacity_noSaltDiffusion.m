    
function c_temp = heatCapacity_noSaltDiffusion(ground)
    
    mineral = ground.STATVAR.mineral;
    porosity = ground.STATVAR.porosity;
    organic = ground.STATVAR.organic;
    a = ground.PARA.a;
    b = ground.PARA.b;
    Tmelt = ground.STATVAR.Tmelt; %in Kelvin
    Tmelt_inDegreeC = Tmelt - 273.15;%convert to °C to compare against temperatures.
    
    T = ground.STATVAR.T;

    %Heat Capacity lives on the midpoints of the cells
    midptDepth = ground.STATVAR.upperPos - ground.STATVAR.layerThick(1)/2 - ...
    cumsum(ground.STATVAR.layerThick);

    %read constants from ground struct
    % J / (K m^3)
    c_water = ground.CONST.c_w;
    c_organic = ground.CONST.c_o;
    c_mineral = ground.CONST.c_m;
    c_ice = ground.CONST.c_i;
    
    L = ground.CONST.L_f;
    
    %determine bulk conductivity and capacity
    c_temp = zeros(size(ground.STATVAR.layerThick));

    %assign heat capacity
    for i = 1:length(midptDepth)
        a1 = 1/porosity(i) - a(i)* abs(Tmelt(i))^b(i);

        %determine water content
        if T(i) <= Tmelt_inDegreeC(i)
            liqWater = 1/(a1 + a(i)*abs(T(i))^b(i));
        else
            liqWater = porosity(i);
        end

        if T(i) <= Tmelt_inDegreeC(i)
            d_liqWater = a(i) * b(i) * abs(T(i)^(b(i)-1) / ...
                (a1 + a(i)*abs(T(i))^b(i))) / (a1 + a(i)*abs(T(i))^b(i));
        else
            d_liqWater = 0;
        end  

        c_temp(i) = c_mineral * mineral(i) + c_water*liqWater + ...
            c_ice*(porosity(i)-liqWater) + c_organic*organic(i) +  L*d_liqWater;
    end  