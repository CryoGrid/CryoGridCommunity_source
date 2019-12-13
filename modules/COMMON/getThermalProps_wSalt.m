function ground = getThermalProps_wSalt(ground)
   
    %read constants
    c_w = ground.CONST.c_w; 
    c_o = ground.CONST.c_o;
    c_m = ground.CONST.c_m;
    c_i = ground.CONST.c_i;
    
    L_f = ground.CONST.L_f;
    tau = ground.CONST.tau;
    
    k_a = ground.CONST.k_a;
    k_w = ground.CONST.k_w;
    k_i = ground.CONST.k_i;
    k_o = ground.CONST.k_o;
    k_m = ground.CONST.k_m;

    
    
    layerThick = ground.STATVAR.layerThick;
    
    
    
    saltDiff0 = ground.PARA.saltDiff0;

    %Calculate liqWater and derivatives
    [liqWater, Tmelt] = watercontent_temperature_salt(ground); %Tmelt in Kelvin!
    ground.STATVAR.liqWater = liqWater;
    ground.STATVAR.Tmelt = Tmelt;
    
    porosity = ground.STATVAR.water;
    
    %ice = ground.STATVAR.ice;
    ice = porosity - liqWater;
        
    mineral = ground.STATVAR.mineral;
    organic = ground.STATVAR.organic;
    air = 1 -liqWater-ice-mineral-organic;
    
    
    [dliqWater_dT, ~] = getDerivativeWaterContent(ground);

    c = c_m.*mineral + c_o.*organic + c_w.*liqWater + c_i .*ice;
    c_eff = c + L_f .* dliqWater_dT;
    %beta=(1 + ns./liqWater.*c./c_eff.*dliqWater_dn);
    %gamma=c_eff + ns./liqWater.*c.*dliqWater_dn;

    saltDiff = saltDiff0.*liqWater./tau;
%     saltDiff = [saltDiff(1);...
%         (saltDiff(1:end-1).*saltDiff(2:end).*(layerThick(1:end-1)+layerThick(2:end)))...
%         ./(saltDiff(1:end-1).*layerThick(2:end)+saltDiff(2:end).*layerThick(1:end-1)); ...
%         saltDiff(end)];
    saltDiff=[saltDiff(1); (saltDiff(1:end-1).*saltDiff(2:end))./(saltDiff(1:end-1)+saltDiff(2:end)); saltDiff(end)];


    thermCond = (liqWater.* k_w.^0.5 + ice.* k_i.^0.5 + mineral.* k_m.^0.5 + organic.* k_o.^0.5 + air.* k_a.^0.5).^2;
    %K=[K(1); (K(1:end-1).*K(2:end))./(K(1:end-1)+K(2:end)); K(end)];
    %The above line is the original code before the K_delta and cT_delta errors
    %were found. This line was missing a factor of 2 which offset the
    %cT_delta/2 error.
    %K=[K(1);((K(1:end-1).*K(2:end)).*(midptDepth(1:end-1)+midptDepth(2:end)))./(K(1:end-1).*midptDepth(2:end)+K(2:end).*midptDepth(1:end-1)); K(end)];
    %Found out that I used midptDepth instead of layerDepth
    %di = layerDepth(2:end)-layerDepth(1:end-1);
    thermCond = [thermCond(1);(thermCond(1:end-1).*thermCond(2:end).*(layerThick(1:end-1)+layerThick(2:end)))./(thermCond(1:end-1).*layerThick(2:end)+thermCond(2:end).*layerThick(1:end-1)); thermCond(end)];
    %K=[K(1);((K(1:end-1).*K(2:end)).*(layerDepth(2:end-1)+layerDepth(3:end)))./(K(1:end-1).*layerDepth(3:end)+K(2:end).*layerDepth(2:end-1)); K(end)];


  
    %update struct
    ground.STATVAR.saltDiff = saltDiff;
    ground.STATVAR.thermCond = thermCond;
    ground.STATVAR.c_eff = c_eff;
end
