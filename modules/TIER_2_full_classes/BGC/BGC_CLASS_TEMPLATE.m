%========================================================================
% CryoGrid GROUND class BGC_CLASS_TEMPLATE

% S. Westermann, November 2020
%========================================================================

classdef BGC_CLASS_TEMPLATE < INITIALIZE

    properties
        PARENT
        IA_BGC
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
            ground.STATVAR.layerThick = [0.01; 0.2; 0.2];
            ground.STATVAR.upperPosition_gridCell = [0; 0.01; 0.21]; %upper position of each grid cell, not used now, but will get relevant for Xice
            
            %initialize with some random values, is overwritten in first timestep
            ia_BGC.BGC.STATVAR.vol_water = ground.STATVAR.layerThick .*0 + 0.5;
            ia_BGC.BGC.STATVAR.vol_mineral = ground.STATVAR.layerThick .*0 + 0.3;
            ia_BGC.BGC.STATVAR.porosity = ground.STATVAR.layerThick .*0 + 0.5;
            ia_BGC.BGC.STATVAR.field_capacity = ground.STATVAR.layerThick .*0 + 0.3;
            ia_BGC.BGC.STATVAR.T = ground.STATVAR.layerThick .*0 + 5;
        end
        
        function ground = provide_CONST(ground)
            %initialize the CONST here
        end
        
        function ground = finalize_init(ground, forcing) 
            %do everything else that is needed to get 
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            %only needs to be done if BGC is doing anything in this timestep
            get_ground_variables(ground.IA_BGC, tile); %get the physical variables over to the BGC class, so that it starts with correct state from the last timestep
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
            
            send_BGC_variables(ground.IA_BGC, tile); %send the BGC variables over to the 
        end
        
        function ground = check_trigger(ground, forcing)
            %add here
        end
        

    end
    
end
