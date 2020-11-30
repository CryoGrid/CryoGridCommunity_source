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
        
        function base = initialize_excel(base)
        
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
        
        function ground = lateral3D_pull_heat(ground, lateral)

        end
        
        function ground = lateral3D_push_heat(ground, lateral)
            
        end
        

    end
    
end

