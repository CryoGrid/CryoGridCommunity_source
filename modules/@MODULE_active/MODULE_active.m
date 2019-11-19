% inhertits from GROUND_base_class and adds SEB as upper boundary
% no interaction with snow is possible here

classdef MODULE_active < GROUND_base_class

     properties
        PARA
     end
    
    
    methods
        
        %mandatory functions for each class

        function ground = provide_variables(ground)  %initializes the subvariables as empty arrays
            
            ground = provide_variables@GROUND_base_class(ground); %call function of the base class
            ground = provide_PARA(ground); %add additional variables
        end
        
%         function ground = assign_global_variables(ground, forcing)
%             ground = assign_global_variables@GROUND_base_class(ground, forcing); %call function of the base class
%             ground.PARA.airT_height = forcing.PARA.airT_height;
%         end
%         
%         function ground = initialize_STATVAR_from_file(ground, grid, forcing, depths)
%             ground = initialize_STATVAR_from_file@GROUND_base_class(ground, grid, forcing, depths);
% 
%             ground = finalize_STATVAR(ground); %assign all variables, that must be calculated or assigned otherwise
%         end
        
    
    end
    
end
