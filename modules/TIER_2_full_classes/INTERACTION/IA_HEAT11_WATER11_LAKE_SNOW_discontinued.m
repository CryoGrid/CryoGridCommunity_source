classdef IA_HEAT11_WATER11_LAKE_SNOW < IA_WATER & IA_HEAT 
    %lake on top of snow - not completed!!
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_m(ia_heat_water);

            get_boundary_condition_BUCKET_LAKE_SNOW_m(ia_heat_water);
        end
        

        
        function trigger_create_LAKE(ia_heat_water, ground, forcing)
            
        end
        

        
        function trigger_remove_LAKE(ia_heat_water, forcing)
            lake = ia_heat_water.PREVIOUS;
            snow = ia_heat_water.NEXT;
            
            snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + sum(lake.STATVAR.waterIce,1);
            snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + sum(lake.STATVAR.energy,1);
            snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) +  sum(lake.STATVAR.waterIce ./ lake.STATVAR.area ,1);
            if sum(strcmp('CHILD', fieldnames(lake))) %unusual case, lake gets removed when it has a snow child - not sure this works
                if lake.CHILD ~= 0
                    snow.STATVAR.waterIce(1) = snow.STATVAR.waterIce(1) + lake.CHILD.STATVAR.waterIce;
                    snow.STATVAR.energy(1) = snow.STATVAR.energy(1) + lake.CHILD.STATVAR.energy;
                    snow.STATVAR.layerThick(1) = snow.STATVAR.layerThick(1) +  lake.CHILD.STATVAR.waterIce ./ lake.CHILD.STATVAR.area;
                end
            end
            
            snow.PREVIOUS = lake.PREVIOUS;
            snow.PREVIOUS.NEXT = snow;
            if ~strcmp(class(snow.PREVIOUS), 'Top')
                snow.IA_PREVIOUS = get_IA_class(class(snow), class(snow.PREVIOUS));
                snow.PREVIOUS.IA_NEXT = snow.IA_PREVIOUS;
                snow.IA_PREVIOUS.PREVIOUS = snow.PREVIOUS;
                snow.IA_PREVIOUS.NEXT = snow;
            end
           snow = compute_diagnostic(snow,forcing);
        end

    end
end