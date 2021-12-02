%========================================================================
% CryoGrid INTERACTION (IA) class for heat conduction and water fluxes between a LAKE class
% and a GROUND class with bucket water and excess ice
% S. Westermann, October 2020
%========================================================================

classdef IA_HEAT11_WATER11_LAKE_XICE < IA_WATER & IA_HEAT
    
    methods
        
        function get_boundary_condition_m(ia_heat_water, tile)
            get_boundary_condition_HEAT_LAKE_m(ia_heat_water);
            get_boundary_condition_BUCKET_LAKE_XICE_m(ia_heat_water); %gravity-driven downwards flow
            get_boundary_condition_BUCKET_LAKE_XWATER_UP_m(ia_heat_water); %Xwater upward flow
        end
        
        %trigger function creating the LAKE class, called by the GROUND
        %class
        function trigger_create_LAKE(ia_heat_water, ground, tile)
                        
%             CURRENT = ground.PREVIOUS;  %go to Top() and get the stored SLEEPING classes
%             while ~strcmp(class(CURRENT), 'Top')
%                 CURRENT = CURRENT.PREVIOUS;
%             end
%             for i=1:size(CURRENT.STORE.SLEEPING,1)  %find the correct sleeping class in the list and copy
%                 if strcmp(class(CURRENT.STORE.SLEEPING{i,1}), ground.PARA.threshold_Xwater_class) && CURRENT.STORE.SLEEPING{i,2} == ground.PARA.threshold_Xwater_index
%                     new_lake = copy(CURRENT.STORE.SLEEPING{i,1});
%                 end
%             end

            for i=1:size(tile.STORE.SLEEPING,1)  %find the correct sleeping class in the list and copy
                if strcmp(class(tile.STORE.SLEEPING{i,1}), ground.PARA.threshold_Xwater_class) && tile.STORE.SLEEPING{i,2} == ground.PARA.threshold_Xwater_index
                    new_lake = copy(tile.STORE.SLEEPING{i,1});
                end
            end
            %initialize state variables
            new_lake.STATVAR.waterIce = ground.STATVAR.Xwater(1);
            new_lake.STATVAR.energy = ground.STATVAR.Xwater(1) .* ground.CONST.c_w .* ground.STATVAR.T(1);
            new_lake.STATVAR.mineral = 0;
            new_lake.STATVAR.organic = 0;
            new_lake.STATVAR.Lstar = ground.STATVAR.Lstar;
            new_lake.STATVAR.area = ground.STATVAR.area(1);
            new_lake.STATVAR.top_depth_rel2groundSurface = 0;
            new_lake.PARA.threshold_Xwater = ground.PARA.threshold_Xwater;
            new_lake.PARA.target_grid = ground.PARA.target_grid;
            new_lake.PARA.target_layerThick = ground.PARA.target_grid;
            
            %subtract the moved part of the state variables
            ground.STATVAR.XwaterIce(1) = ground.STATVAR.XwaterIce(1) - new_lake.STATVAR.waterIce;
            ground.STATVAR.Xwater(1) = ground.STATVAR.Xwater(1) - new_lake.STATVAR.waterIce;
            ground.STATVAR.energy(1) = ground.STATVAR.energy(1) - new_lake.STATVAR.energy;
            ground.STATVAR.layerThick(1) = ground.STATVAR.layerThick(1) -  new_lake.STATVAR.waterIce ./ new_lake.STATVAR.area;
            
            % if CHILD snow exists, add the CHILD snow to the new LAKE and remove CHILD
            if sum(strcmp('CHILD', fieldnames(ground)))
                if ground.CHILD ~= 0
                    new_lake.STATVAR.waterIce = new_lake.STATVAR.waterIce + ground.CHILD.STATVAR.waterIce;
                    new_lake.STATVAR.energy = new_lake.STATVAR.energy + ground.CHILD.STATVAR.energy;
                    ground.CHILD = 0;
                    ground.IA_CHILD = 0;
                end
            end
            new_lake.STATVAR.layerThick = new_lake.STATVAR.waterIce ./ new_lake.STATVAR.area;
            new_lake = compute_diagnostic(new_lake, tile);
            
            %change stratigraphy
            ground.PREVIOUS.NEXT = new_lake;
            new_lake.NEXT = ground;
            new_lake.PREVIOUS = ground.PREVIOUS;
            ground.PREVIOUS = new_lake;
            
            %ground.PREVIOUS.NEXT.NEXT = ground;
            ground.IA_PREVIOUS = get_IA_class(class(new_lake), class(ground)); %dedicated function which calls the TIER2 function
            ground.PREVIOUS.IA_NEXT = ground.IA_PREVIOUS;
            ground.IA_PREVIOUS.PREVIOUS = new_lake;
            ground.IA_PREVIOUS.NEXT = ground;
        end
        
                
        %trigger function removing the LAKE class, called by the LAKE 
        %class itself
        function trigger_remove_LAKE(ia_heat_water, tile)
            lake = ia_heat_water.PREVIOUS;
            ground = ia_heat_water.NEXT;
            
            ground.STATVAR.XwaterIce(1) = ground.STATVAR.XwaterIce(1) + sum(lake.STATVAR.waterIce,1);
            
            %changed Sep 2020
            E_frozen = -sum(lake.STATVAR.waterIce,1).* ground.CONST.L_f;
            unfrozen_fraction = max(0, min(1, (sum(lake.STATVAR.energy,1)-E_frozen)./-E_frozen));
            ground.STATVAR.Xwater(1) = ground.STATVAR.Xwater(1) + sum(lake.STATVAR.waterIce,1) .* unfrozen_fraction; %sum(lake.STATVAR.waterIce,1);
            
            
            ground.STATVAR.energy(1) = ground.STATVAR.energy(1) + sum(lake.STATVAR.energy,1);
            ground.STATVAR.layerThick(1) = ground.STATVAR.layerThick(1) +  sum(lake.STATVAR.waterIce ./ lake.STATVAR.area ,1);
            if sum(strcmp('CHILD', fieldnames(lake))) %pass snow child down from lake to ground
                if lake.CHILD ~= 0
                    ground.CHILD = lake.CHILD;
                    ground.CHILD.PARENT = ground;
                    ground.CHILD.NEXT = ground; 
                    ground.IA_CHILD = get_IA_class(class(ground.CHILD), class(ground));
                    ground.IA_CHILD.PREVIOUS = ground.CHILD;
                    ground.IA_CHILD.NEXT = ground; %SNOW CHILD created
                end
            end
            
            ground.PREVIOUS = lake.PREVIOUS;
            ground.PREVIOUS.NEXT = ground;
            if ~strcmp(class(ground.PREVIOUS), 'Top')
                ground.IA_PREVIOUS = get_IA_class(class(ground.PREVIOUS), class(ground));
                ground.PREVIOUS.IA_NEXT = ground.IA_PREVIOUS;
                ground.IA_PREVIOUS.PREVIOUS = ground.PREVIOUS;
                ground.IA_PREVIOUS.NEXT = ground;
            end
           ground = compute_diagnostic(ground,tile);
        end
                
    end
end