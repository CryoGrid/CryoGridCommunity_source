function ground = get_boundary_salt_u(ground, forcing)

if forcing.TEMP.saltConcForcing == 0
    saltFlux_ub = 0;
else
    saltFlux_ub = ground.STATVAR.saltDiff(1).*(ground.STATVAR.saltConc(1) - forcing.TEMP.saltConcForcing)./ ...
        abs(ground.STATVAR.layerThick(1)/2);%abs(ground.STATVAR.midptDepth(1) - ground.STATVAR.layerDepth(1));
end



ground.TEMP.saltFlux_ub = saltFlux_ub;
ground.TEMP.saltConc_ub = forcing.TEMP.saltConcForcing;

