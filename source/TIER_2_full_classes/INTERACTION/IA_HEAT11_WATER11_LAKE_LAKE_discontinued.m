classdef IA_HEAT11_WATER11_LAKE_LAKE < IA_WATER & IA_HEAT 
    
    methods
        
        function get_boundary_condition_m(ia_heat_water)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_LAKE_LAKE_m(ia_heat_water);
        end
        

        
        function trigger_remove_LAKE(ia_heat_water, forcing)
            lake1 = ia_heat_water.PREVIOUS;
            lake2 = ia_heat_water.NEXT;
            
            lake2.STATVAR.waterIce(1) = lake2.STATVAR.XwaterIce(1) + sum(lake1.STATVAR.waterIce,1);
            lake2.STATVAR.energy(1) = lake2.STATVAR.energy(1) + sum(lake1.STATVAR.energy,1);
            lake2.STATVAR.layerThick(1) = lake2.STATVAR.layerThick(1) +  sum(lake1.STATVAR.waterIce ./ lake1.STATVAR.area ,1);
            if sum(strcmp('CHILD', fieldnames(lake1)))
                if lake2.CHILD ~= 0
                    lake2.CHILD = lake1.CHILD;
                    lake2.CHILD.PARENT = lake2;
                    lake2.CHILD.NEXT = lake2; 
                    lake2.IA_CHILD = get_IA_class(class(lake2.CHILD), class(lake2));
                    lake2.IA_CHILD.PREVIOUS = lake2.CHILD;
                    lake2.IA_CHILD.NEXT = lake2; %SNOW CHILD created
                end
            end
            
            lake2.PREVIOUS = lake1.PREVIOUS;
            lake2.PREVIOUS.NEXT = lake2;
            if ~strcmp(class(lake2.PREVIOUS), 'Top')
                lake2.IA_PREVIOUS = get_IA_class(class(lake2), class(lake2.PREVIOUS));
                lake2.PREVIOUS.IA_NEXT = lake2.IA_PREVIOUS;
                lake2.IA_PREVIOUS.PREVIOUS = lake2.PREVIOUS;
                lake2.IA_PREVIOUS.NEXT = lake2;
            end
           lake2 = compute_diagnostic(lake2,forcing);
        end
            
           

    end
end