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
                
                lateral.PARENT.STATVAR2ALL.ds_top_snow_date = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.top_snow_date) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_top_snow_date(isnan(lateral.PARENT.STATVAR2ALL.ds_top_snow_date)) = 0;
                lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.bottom_snow_date) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date(isnan(lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date)) = 0;
            else
                lateral.PARENT.STATVAR2ALL.snow_drift = 1; 
            end
        end
        
        function snow = lateral3D_pull_snow_crocus_dump(snow, lateral)
            snow = prog_wind_drift(snow); %populate one_over_tau
            
            %fraction_mobile = snow.TEMP.one_over_tau .* lateral.PARA.N_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment;
            
            drifting_SWE =  snow.TEMP.one_over_tau .* lateral.PARA.n_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment .* snow.STATVAR.area .* lateral.PARA.weighting_factor; %absolute rate, potentially driftable snow
            %-lateral.STATVAR.exposure .* lateral.STATVAR.ratio_receiving_total_area
            fraction_mobile = drifting_SWE ./ snow.STATVAR.ice;
            fraction_mobile(isnan(fraction_mobile))=0;
            fraction_mobile = min(0.5, fraction_mobile);
            
            lateral.PARENT.STATVAR2ALL.upperPos = snow.STATVAR.upperPos;
            lateral.PARENT.STATVAR2ALL.ds_area = snow.STATVAR.area(1,1);
            lateral.PARENT.STATVAR2ALL.weighting_factor = lateral.PARA.weighting_factor;
            if max(fraction_mobile) > 0 && sum(snow.STATVAR.layerThick,1) > lateral.PARA.snow_holding_height 
                lateral.PARENT.STATVAR2ALL.snow_drift = 2;
                %fraction_mobile = min(0.5, fraction_mobile);
                
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
                
                lateral.PARENT.STATVAR2ALL.ds_top_snow_date = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.top_snow_date) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_top_snow_date(isnan(lateral.PARENT.STATVAR2ALL.ds_top_snow_date)) = 0;
                lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date = sum(fraction_mobile .* snow.STATVAR.waterIce .* snow.STATVAR.bottom_snow_date) ./ lateral.PARENT.STATVAR2ALL.ds_waterIce;
                lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date(isnan(lateral.PARENT.STATVAR2ALL.ds_bottom_snow_date)) = 0;
            else
                lateral.PARENT.STATVAR2ALL.snow_drift = 1; 
            end
        end
        
        %push function for SNOW_crocus_...
        function snow = lateral3D_push_snow_crocus(snow, lateral)
            if lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure < 0  && lateral.PARENT.STATVAR2ALL.snow_drift == 2 %loose snow
                remaining_fraction = 1 - min(0.5, snow.TEMP.one_over_tau .* lateral.PARA.N_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment);

                snow.STATVAR.energy = remaining_fraction.* snow.STATVAR.energy;
                snow.STATVAR.waterIce = remaining_fraction.* snow.STATVAR.waterIce;
                snow.STATVAR.water = remaining_fraction.* snow.STATVAR.water;
                snow.STATVAR.ice = remaining_fraction.* snow.STATVAR.ice;
                snow.STATVAR.layerThick = remaining_fraction.* snow.STATVAR.layerThick;
                %in principle to_snow_date should be changed as well, but
                %probably not critical
            elseif lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure > 0 %gain snow
                new_snow.STATVAR = lateral.STATVAR.ds;
                
                snow = merge_cells_intensive2(snow, 1, new_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'; 'target_density';}, 'ice');
                snow = merge_cells_extensive2(snow, 1, new_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'; 'water'});
                snow = merge_cells_snowfall_times2(snow, 1, new_snow, 1); %specific function merginging bottom and top snow dates
            end
        end
        
        %push function for SNOW_crocus2_...
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
                snow = merge_cells_snowfall_times2(snow, 1, new_snow, 1); %specific function merginging bottom and top snow dates
                
                snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1);
                snow.STATVAR.target_density(1) = snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1);
                snow.STATVAR.target_density(1) = min(1, snow.STATVAR.target_density(1));  % avoids rounding errors, if >1, this might trigger a follow-up problem
                snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1) ./ snow.STATVAR.area(1,1));
            end
        end
        
        function snow = lateral3D_push_snow_crocus_dump(snow, lateral)
            if lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure < 0  && lateral.PARENT.STATVAR2ALL.snow_drift == 2 %loose snow
                drifting_SWE = -lateral.STATVAR.exposure .* lateral.STATVAR.ratio_receiving_total_area .* snow.TEMP.one_over_tau .* lateral.PARA.n_drift .* ...
                    lateral.CONST.day_sec .* lateral.PARA.ia_time_increment .* snow.STATVAR.area .* lateral.PARA.weighting_factor; %absolute rate
                drifting_SWE = min(drifting_SWE, 0.5 .* snow.STATVAR.ice); %limit amount of drifting snow per timestep for numerical reasons
                remaining_fraction = 1 - drifting_SWE ./ snow.STATVAR.ice;

                snow.STATVAR.energy = remaining_fraction.* snow.STATVAR.energy;
                snow.STATVAR.waterIce = remaining_fraction.* snow.STATVAR.waterIce;
                snow.STATVAR.water = remaining_fraction.* snow.STATVAR.water;
                snow.STATVAR.ice = remaining_fraction.* snow.STATVAR.ice;
                snow.STATVAR.layerThick = remaining_fraction.* snow.STATVAR.layerThick;

            elseif lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure > 0 %gain snow
                new_snow.STATVAR = lateral.STATVAR.ds;
                
                snow = merge_cells_intensive2(snow, 1, new_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'; 'target_density';}, 'ice');
                snow = merge_cells_extensive2(snow, 1, new_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'; 'water'});
                snow = merge_cells_snowfall_times2(snow, 1, new_snow, 1); %specific function merginging bottom and top snow dates
            end
        end
        
        %push function for SNOW_crocus2_...
        function snow = lateral3D_push_snow_crocus2_dump(snow, lateral)
            
            if lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure < 0  && lateral.PARENT.STATVAR2ALL.snow_drift == 2 %loose snow
%                 drifting_SWE = -lateral.STATVAR.exposure .* lateral.STATVAR.ratio_receiving_total_area .* snow.TEMP.one_over_tau .* lateral.PARA.n_drift .* ...
%                     lateral.CONST.day_sec .* lateral.PARA.ia_time_increment .* snow.STATVAR.area .* lateral.PARA.weighting_factor; %absolute rate
                
                drifting_SWE =  snow.TEMP.one_over_tau .* lateral.PARA.n_drift .* lateral.CONST.day_sec .* lateral.PARA.ia_time_increment .* snow.STATVAR.area .* lateral.PARA.weighting_factor; %absolute rate, potentially driftable snow
                %-lateral.STATVAR.exposure .* lateral.STATVAR.ratio_receiving_total_area
                fraction_mobile = drifting_SWE ./ snow.STATVAR.ice;
                fraction_mobile(isnan(fraction_mobile))=0;
                fraction_mobile = min(0.5, fraction_mobile);
                drifting_SWE = -lateral.STATVAR.exposure .* fraction_mobile .* snow.STATVAR.ice;
                lateral.STATVAR.drift_pool = lateral.STATVAR.drift_pool - sum(drifting_SWE);
                
                %drifting_SWE = min(drifting_SWE, 0.5 .* snow.STATVAR.ice); %limit amount of drifting snow per timestep for numerical reasons
                remaining_fraction = 1 - drifting_SWE ./ snow.STATVAR.ice;
                
                snow.STATVAR.energy = remaining_fraction.* snow.STATVAR.energy;
                snow.STATVAR.waterIce = remaining_fraction.* snow.STATVAR.waterIce;
                snow.STATVAR.water = remaining_fraction.* snow.STATVAR.water;
                snow.STATVAR.ice = remaining_fraction.* snow.STATVAR.ice;
                snow.STATVAR.layerThick = remaining_fraction.* snow.STATVAR.layerThick;
                
            elseif lateral.STATVAR.snow_drift_yes_no && lateral.STATVAR.exposure > 0 %gain snow
                new_snow.STATVAR = lateral.STATVAR.ds;
                lateral.STATVAR.drift_pool = lateral.STATVAR.drift_pool + lateral.STATVAR.ds.waterIce;
                snow.STATVAR.layerThick(1) =  snow.STATVAR.layerThickSnowFirstCell;
                
                snow = merge_cells_intensive2(snow, 1, new_snow, 1, {'d'; 's'; 'gs'; 'time_snowfall'}, 'ice');
                snow = merge_cells_extensive2(snow, 1, new_snow, 1, {'waterIce'; 'energy'; 'layerThick'; 'ice'});
                snow = merge_cells_snowfall_times2(snow, 1, new_snow, 1); %specific function merginging bottom and top snow dates
                
                snow.STATVAR.layerThickSnowFirstCell = snow.STATVAR.layerThick(1);
                snow.STATVAR.target_density(1) = snow.STATVAR.ice(1) ./ snow.STATVAR.layerThick(1) ./ snow.STATVAR.area(1);
                snow.STATVAR.target_density(1) = min(1, snow.STATVAR.target_density(1));  % avoids rounding errors, if >1, this might trigger a follow-up problem
                snow.STATVAR.layerThick(1) = max(snow.STATVAR.layerThick(1), snow.STATVAR.waterIce(1) ./ snow.STATVAR.area(1,1));
            end
        end
        

    end
end

