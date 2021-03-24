%========================================================================
% CryoGrid INTERACTION (IA) class connecting the frozen and unfrozen
% classes of LAKE_simple
% NOTE: also contains experimental code for boundary conditions and triggers between two LAKE classes 
% NOTE: used for LAKE classes with and without wtaer cycle
% S. Westermann, October 2020
%========================================================================

classdef IA_LAKE_simple_frozen_unfrozen < IA_WATER & IA_HEAT 
    
    methods

        function lake_next_season = create_annihilate(ia_create_next_season_lake, lake_this_season)
            class_handle = str2func(lake_this_season.PARA.next_season_lake_class);
            %lake_next_season = class_handle(-1,0,0,0); 
            lake_next_season = class_handle(); 
            lake_next_season = initialize_from_LAKE_previous_season(lake_next_season, lake_this_season);
        end
        

        %experimental code connecting two LAKES, not used at this stage
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_m(ia_heat_water);
            get_boundary_condition_BUCKET_LAKE_LAKE_m(ia_heat_water);
        end
        
        
        function trigger_remove_LAKE(ia_heat_water, tile)
            lake1 = ia_heat_water.PREVIOUS;
            lake2 = ia_heat_water.NEXT;
            
            lake2.STATVAR.waterIce(1) = lake2.STATVAR.waterIce(1) + sum(lake1.STATVAR.waterIce,1);
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
           lake2 = compute_diagnostic(lake2,tile);
        end
        
    end
end