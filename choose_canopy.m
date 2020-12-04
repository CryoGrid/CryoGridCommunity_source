function ground = choose_canopy(ground, t)
   
%     if t>735882 && t<736064 || t>736247 && t<736430 || t>736613 && t<736795  || t>736978 && t<737160 || t>737343 && t<737525
    % 30.5 - 10.10 (summer period)
    if t>735882 && t<736114 || t>736247 && t<736480 || t>736613 && t<736845  || t>736978 && t<737210 || t>737343 && t<737575
        ground.STATVAR.vegetation.canopy.pai = ground.STATVAR.vegetation.canopy.pai_winter;
        ground.STATVAR.vegetation.canopy.lai = ground.STATVAR.vegetation.canopy.lai_winter;
        ground.STATVAR.vegetation.canopy.dlai = ground.STATVAR.vegetation.canopy.dlai_winter(:,:);
        ground.STATVAR.vegetation.canopy.sumlai = ground.STATVAR.vegetation.canopy.sumlai_winter(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.zw = ground.STATVAR.vegetation.mlcanopyinst.zw_winter(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.zs = ground.STATVAR.vegetation.mlcanopyinst.zs_winter(:,:);
        ground.STATVAR.vegetation.canopy.dpai = ground.STATVAR.vegetation.canopy.dpai_winter(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.sumpai = ground.STATVAR.vegetation.mlcanopyinst.sumpai_winter(:,:);
        ground.STATVAR.vegetation.flux.albsoib = ground.STATVAR.vegetation.flux.albsoib_winter(:,:);
        ground.STATVAR.vegetation.flux.albsoid = ground.STATVAR.vegetation.flux.albsoid_winter(:,:);
    else % Summer canopy
        ground.STATVAR.vegetation.canopy.pai = ground.STATVAR.vegetation.canopy.pai_summer;
        ground.STATVAR.vegetation.canopy.lai = ground.STATVAR.vegetation.canopy.lai_summer;
        ground.STATVAR.vegetation.canopy.dlai = ground.STATVAR.vegetation.canopy.dlai_summer(:,:);
        ground.STATVAR.vegetation.canopy.sumlai = ground.STATVAR.vegetation.canopy.sumlai_summer(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.zw = ground.STATVAR.vegetation.mlcanopyinst.zw_summer(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.zs = ground.STATVAR.vegetation.mlcanopyinst.zs_summer(:,:);
        ground.STATVAR.vegetation.canopy.dpai = ground.STATVAR.vegetation.canopy.dpai_summer(:,:);
        ground.STATVAR.vegetation.mlcanopyinst.sumpai = ground.STATVAR.vegetation.mlcanopyinst.sumpai_summer(:,:);
        ground.STATVAR.vegetation.flux.albsoib = ground.STATVAR.vegetation.flux.albsoib_summer(:,:);
        ground.STATVAR.vegetation.flux.albsoid = ground.STATVAR.vegetation.flux.albsoid_summer(:,:);
    end
    
    ground.STATVAR.vegetation.pftcon.vcmaxpft = 51;
    ground.STATVAR.vegetation.pftcon.slatop = 0.024; 
      
end