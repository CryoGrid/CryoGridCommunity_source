%========================================================================
% CryoGrid GROUND class BGC_CLASS_TEMPLATE

% S. Westermann, November 2020
%========================================================================

classdef BGC_CLASS_TEMPLATE < INITIALIZE

    properties
        PARENT
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = BGC_CLASS_TEMPLATE(index, pprovider, cprovider, forcing)  
            ground@INITIALIZE(index, pprovider, cprovider, forcing);
        end
        
         function ground = provide_PARA(ground)
            %initialize the PARA here
            
        end
        
        function ground = provide_STATVAR(ground)
            %initialize the STATVAR here
            GROUND.STATVAR.layerThick = [0.2; 0.2];
            ground.STATVAR.upperPosition_gridCell = [0; 0.2]; %upper position of each grid cell, normally [0; cumsum(GROUND.STATVAR.layerThick(1:end-1))]
        end
        
        function ground = provide_CONST(ground)
            %initialize the CONST here
        end
        
        function ground = finalize_init(ground, forcing) 
            %do everything else that is needed to get 
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            %add here
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
             %add here
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            %add here
        end
        
        function timestep = get_timestep(ground, tile)
            %modify
            timestep = 1e9;
        end
        
        function ground = advance_prognostic(ground, tile) 
            %add here
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            %add here
        end
       
        function ground = compute_diagnostic(ground, tile)
            %add here
        end
        
        function ground = check_trigger(ground, forcing)
            %add here
        end
        

    end
    
end
