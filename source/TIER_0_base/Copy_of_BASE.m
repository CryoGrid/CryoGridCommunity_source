%========================================================================
% CryoGrid BASE class for all GROUND classes 
% S. Westermann, Oct 2020
%========================================================================

classdef BASE < matlab.mixin.Copyable
    
    properties
        class_index
        CONST       %constants
        PARA        %external service parameters
        STATVAR     %state varibales, e.g. energy, water content
        TEMP        %mainly derivatives in prognostic timestep 
        PREVIOUS    %pointer to GROUND class above current class in stratigraphy
        NEXT        %pointer to GROUND class below current class in stratigraphy
        IA_PREVIOUS  %pointer to IA class with GROUND class above current class in stratigraphy
        IA_NEXT      %pointer to IA class with GROUND class below current class in stratigraphy
    end
    
    methods 
        
        %init_steady_state, makes a steady_state T profile
        function [ground, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = init_T_steady_state_TOP_CLASS(ground, T_first_cell, start_depth_steady_state, heat_flux)
            ground = finalize_init2(ground);
            i=1;
            depth = start_depth_steady_state; %find correct start cell
            while i<= size(ground.STATVAR.T,1) && depth - ground.STATVAR.layerThick(i,1) >= 0
                depth = depth - ground.STATVAR.layerThick(i,1);
                i=i+1;
            end
            if i==size(ground.STATVAR.T,1) && depth > 0
                start_depth_steady_state = depth;
            else
                for j=i+1:size(ground.STATVAR.T,1)
                    ground.STATVAR.T(j,1) =  ground.STATVAR.T(j-1,1) + heat_flux ./ ground.STATVAR.thermCond(j-1,1) .* (ground.STATVAR.layerThick(j-1,1) + ground.STATVAR.layerThick(j,1))./2;
                    ground = finalize_init2(ground);
                end
                start_depth_steady_state = NaN;
            end
            
            T_end = ground.STATVAR.T(end,1);
            thermCond_end = ground.STATVAR.thermCond(end,1);
            layerThick_end = ground.STATVAR.layerThick(end,1);
        end
        
        function [ground, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = init_T_steady_state(ground, T_end, start_depth_steady_state, thermCond_end, layerThick_end, heat_flux)
            i=1;
            if isnan(start_depth_steady_state)
                ground.STATVAR.T(1,1) = T_end + heat_flux ./ thermCond_end .* (layerThick_end + ground.STATVAR.layerThick(1,1))./2;
            else
                depth = start_depth_steady_state; %find correct start cell
                while i<= size(ground.STATVAR.T,1) && depth - ground.STATVAR.layerThick(i,1) > 0
                    depth = depth - ground.STATVAR.layerThick(i,1);
                    i=i+1;
                end
            end
            ground = finalize_init2(ground);

            if i==size(ground.STATVAR.T,1) && depth > 0
                start_depth_steady_state = depth;
            else
                for j=i+1:size(ground.STATVAR.T,1)
                    ground.STATVAR.T(j,1) =  ground.STATVAR.T(j-1,1) + heat_flux ./ ground.STATVAR.thermCond(j-1,1) .* (ground.STATVAR.layerThick(j-1,1) + ground.STATVAR.layerThick(j,1))./2;
                    ground = finalize_init2(ground);
                end
                start_depth_steady_state = NaN;
            end
            
            T_end = ground.STATVAR.T(end,1);
            thermCond_end = ground.STATVAR.thermCond(end,1);
            layerThick_end = ground.STATVAR.layerThick(end,1);
        end
        
        %----LATERAL----------
        
        %empty functions for lateral interactions, overwritten in TIER 2-3 if the functions are active
        %all functions used in lateral IA classes must be declared here
        function base = lateral_push_remove_surfaceWater(base, lateral)
            
        end
        
        function base = lateral_push_remove_subsurfaceWater(base, lateral)
            
        end
        
        function ground = lateral_push_remove_water_seepage(ground, lateral)
            
        end
        
        function ground = lateral_push_water_reservoir(ground, lateral)
            lateral.TEMP.open_system = 0;
        end
        
        function ground = lateral_push_remove_water_overland_flow(ground, lateral)
            
        end
        
        function ground = lateral_push_heat(ground, lateral)
            
        end
        
        function ground = lateral3D_pull_water_unconfined_aquifer(ground, lateral)
            lateral.TEMP.open_system = 0;
        end
        
        function ground = lateral3D_push_water_unconfined_aquifer(ground, lateral)
            
        end
        
        function ground = lateral3D_pull_water_general_aquifer(ground, lateral)
            lateral.TEMP.open_system = 0;
        end
        
        function ground = lateral3D_push_water_general_aquifer(ground, lateral)
            
        end
        
        function [saturated_next, hardBottom_next] = get_saturated_hardBottom_first_cell(ground, lateral)
            saturated_next = 0;
            hardBottom_next = 1;
        end
        
        function ground = lateral3D_pull_water_overland_flow(ground, lateral)
            
        end
        
        function ground = lateral3D_pull_heat(ground, lateral)

        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            
        end
        
        function ground = reset_time_BGC(ground, tile)
            
        end
        
        

    end
    
end

