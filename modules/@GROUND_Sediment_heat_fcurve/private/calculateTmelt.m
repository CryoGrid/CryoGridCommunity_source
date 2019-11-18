function [Tmelt, a, b] = calculateTmelt(ground)
    %calculate Tmelt in K
    
    load('lookup_table_fc.mat')
    a_sand=interp1(fc_parameters(:,1), fc_parameters(:,2), ground.STATVAR.saltConc).* ...
        ground.STATVAR.porosity.^(-1);
    b_sand=interp1(fc_parameters(:,1), fc_parameters(:,3), ground.STATVAR.saltConc);
    a_silt=interp1(fc_parameters(:,1), fc_parameters(:,4), ground.STATVAR.saltConc).* ...
        ground.STATVAR.porosity.^(-1);
    b_silt=interp1(fc_parameters(:,1), fc_parameters(:,5), ground.STATVAR.saltConc);
    a = a_silt.*(ground.STATVAR.soilType==1)+a_sand.*(ground.STATVAR.soilType==0);
    b = b_silt.*(ground.STATVAR.soilType==1)+b_sand.*(ground.STATVAR.soilType==0);

    
    Tmelt = ground.CONST.Tmelt_free_water./ground.CONST.L_f.* ...
        (-ground.CONST.R .* ground.STATVAR.saltConc.*ground.CONST.Tmelt_free_water);

    
end

