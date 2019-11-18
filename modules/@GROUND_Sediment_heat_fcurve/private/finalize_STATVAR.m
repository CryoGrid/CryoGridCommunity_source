function ground = finalize_STATVAR(ground)

    midptDepth = ground.STATVAR.upperPos - ground.STATVAR.layerThick(1)/2 - ...
    cumsum(ground.STATVAR.layerThick);

    %calculate porosity
    porosityZero=1-(ground.STATVAR.mineral+ground.STATVAR.organic);
    bulkDensityZero=1./((porosityZero+0.6845)/1.8);
    %bulkDensity =bulkDensityZero-5.1.*1e-5.*D + 0.0037 .*z.^0.766;
    bulkDensity =bulkDensityZero + 0.0037.*abs(midptDepth).^0.766; %abs to avoid imaginary numbers!
    
    porosity=1.80.*bulkDensity.^(-1)-0.6845;
    porosity(porosity<0.03)=0.03;
    
    
    %update statvar
    ground.STATVAR.porosity = porosity;
    ground.STATVAR.mineral=ground.STATVAR.mineral.*(1-porosity)./(ground.STATVAR.mineral+ground.STATVAR.organic);
    ground.STATVAR.organic=ground.STATVAR.organic.*(1-porosity)./(ground.STATVAR.mineral+ground.STATVAR.organic);
    ground.STATVAR.water=1-ground.STATVAR.organic-ground.STATVAR.mineral;

    %calculate state variables that depend on the more constant ones
    [Tmelt, a, b] = calculateTmelt(ground);
    ground.STATVAR.Tmelt = Tmelt; % in Kelvin!
    ground.PARA.a = a;
    ground.PARA.b = b;

    %Calculate initial condition T0
    T = steadyState(ground); %steady state
    %or spinup - put that here
    ground.STATVAR.T = T;  % [degree C]


    %conductivity, heat capacity and liquid water content
    ground = getThermalProps_noSaltDiffusion(ground);
end
