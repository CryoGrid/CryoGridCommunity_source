%========================================================================
% CryoGrid TIER1 library class for functions for lateral fluxes of snow
% only active for CROCUS snow classes at this time
% contains push and pull functions use for LAT3D_SNOW 
% S. Westermann, October 2020
%========================================================================


classdef SNOW_FLUXES_LATERAL < BASE
    
    methods
                
        
        function snow = lateral3D_pull_snow_crocus(snow, lateral)
            snow = prog_wind_drift(snow); %populate one_over_tau
            
            fraction_mobile = snow.TEMP.one_over_tau .* lateral.PARA.N_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment;
            
            lateral.PARENT.STATVAR2ALL.upperPos = snow.STATVAR.upperPos;
            lateral.PARENT.STATVAR2ALL.ds_area = snow.STATVAR.area(1,1);
            if max(fraction_mobile) > 0
                lateral.PARENT.STATVAR2ALL.snow_drift = 2;
                fraction_mobile = min(0.5, fraction_mobile);
                
                %extensive variables
                lateral.PARENT.STATVAR2ALL.ds_volume = sum(fraction_mobile .* snow.STATVAR.layerThick .* snow.STATVAR.area);
                lateral.PARENT.STATVAR2ALL.ds_waterIce = sum(fraction_mobile .* snow.STATVAR.waterIce);
                lateral.PARENT.STATVAR2ALL.ds_energy = sum(fraction_mobile .* snow.STATVAR.energy);
                lateral.PARENT.STATVAR2ALL.ds_ice = sum(fraction_mobile .* snow.STATVAR.ice);
                %intensive variables - use waterIce as scaling variable,
                %identical to ice when snow is driftable
                lateral.PARENT.STATVAR2ALL.ds_d = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.d) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_d(isnan(lateral.PARENT.STATVAR2ALL.ds_d)) = 0;
                lateral.PARENT.STATVAR2ALL.ds_s = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.s) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_s(isnan(lateral.PARENT.STATVAR2ALL.ds_s)) = 0;
                lateral.PARENT.STATVAR2ALL.ds_gs = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.gs) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_gs(isnan(lateral.PARENT.STATVAR2ALL.ds_gs)) = 0;
                lateral.PARENT.STATVAR2ALL.ds_time_snowfall = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.time_snowfall) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_time_snowfall(isnan(lateral.PARENT.STATVAR2ALL.ds_time_snowfall)) = 0;
            else
                lateral.PARENT.STATVAR2ALL.snow_drift = 1; 
            end
        end
        
        function snow = lateral3D_push_snow_crocus(snow, lateral)
            if lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure < 0  && lateral.PARENT.STATVAR2ALL.snow_drift == 2 %loose snow
                remaining_fraction = 1 - min(0.5, snow.TEMP.one_over_tau .* lateral.PARA.N_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment);

                snow.STATVAR.energy = remaining_fraction.* snow.STATVAR.energy;
                snow.STATVAR.waterIce = remaining_fraction.* snow.STATVAR.waterIce;
                snow.STATVAR.water = remaining_fraction.* snow.STATVAR.water;
                snow.STATVAR.ice = remaining_fraction.* snow.STATVAR.ice;
                snow.STATVAR.layerThick = remaining_fraction.* snow.STATVAR.layerThick;

            elseif lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure > 0 %gain snow
                new_snow.STATVAR = lateral.STATVAR.ds;
                
                snow = merge_cells_intensive2(snow, 1, new_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'; 'target_density';}, 'ice');
                snow = merge_cells_extensive2(snow, 1, new_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'; 'water'});
            end
        end
        
        function snow = lateral3D_push_snow_crocus2(snow, lateral)
            
            if lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure < 0  && lateral.PARENT.STATVAR2ALL.snow_drift == 2 %loose snow
                remaining_fraction = 1 - min(0.5, snow.TEMP.one_over_tau .* lateral.PARA.N_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment);
                
                snow.STATVAR.energy = remaining_fraction.* snow.STATVAR.energy;
                snow.STATVAR.waterIce = remaining_fraction.* snow.STATVAR.waterIce;
                snow.STATVAR.water = remaining_fraction.* snow.STATVAR.water;
                snow.STATVAR.ice = remaining_fraction.* snow.STATVAR.ice;
                snow.STATVAR.layerThick = remaining_fraction.* snow.STATVAR.layerThick;
                
            elseif lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure > 0 %gain snow
                new_snow.STATVAR = lateral.STATVAR.ds;
                snow.STATVAR.layerThick(1) =  snow.STATVAR.layerThickSnowFirstCell;
                
                snow = merge_cells_intensive2(snow, 1, new_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow = merge_cells_extensive2(snow, 1, new_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'});
                
                snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1);
                snow.STATVAR.target_density(1) = snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1);
                snow.STATVAR.target_density(1) = min(1, snow.STATVAR.target_density(1));  % avoids rounding errors, if >1, this might trigger a follow-up problem
                snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1) ./ snow.STATVAR.area(1,1));
            end
        end

    end
end

