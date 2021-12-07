%========================================================================
% CryoGrid GROUND class READ_STATVAR_FROM_OUT_BGC
% reads variables from OUT files of the OUT class
% OUT_all_lateral_STORE4READ and interacts with BGC class
% S. Westermann, November 2020
%========================================================================

classdef READ_STATVAR_FROM_OUT_BGC < READ_STATVAR_FROM_OUT

    properties
        BGC
        IA_BGC
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
         function ground = provide_PARA(ground)
            
            ground = provide_PARA@READ_STATVAR_FROM_OUT(ground);

            ground.PARA.BGC_CLASS = [];
            
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@READ_STATVAR_FROM_OUT(ground);
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.day_sec = [];
        end
        
        function ground = finalize_init(ground, tile) 
            ground = finalize_init@READ_STATVAR_FROM_OUT(ground, tile);

            class_handle = str2func(ground.RUN_PARA.BGC_CLASS);
            ground.BGC = class_handle(); 
            ground.BGC.PARENT = ground;
            %remove this in the end
            ground.BGC = provide_PARA(ground.BGC);
            ground.BGC = provide_STATVAR(ground.BGC);
            ground.BGC = provide_CONST(ground.BGC);
            ground.BGC = finalize_init(ground.BGC, tile);
            
            ground.IA_BGC = IA_BGC_read_statvar_from_out();
            ground.IA_BGC.BGC = ground.BGC;
            ground.IA_BGC.GROUND = ground;
            ground.BGC.IA_BGC = ground.IA_BGC;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            ground = get_boundary_condition_u@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = get_boundary_condition_u(ground.BGC, tile);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = get_boundary_condition_l(ground.BGC, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivatives_prognostic@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = get_derivatives_prognostic(ground.BGC, tile);
        end
        
        function timestep = get_timestep(ground, tile)
            timestep_ground = get_timestep@READ_STATVAR_FROM_OUT(ground, tile);
            timestep_BGC = get_timestep(ground.BGC, tile);
            timestep = min(timestep_ground, timestep_BGC);
        end
        
        function ground = advance_prognostic(ground, tile) 
            ground = advance_prognostic@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = advance_prognostic(ground.BGC, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            ground = compute_diagnostic_first_cell@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = compute_diagnostic_first_cell(ground.BGC, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)

            ground.BGC = compute_diagnostic(ground.BGC, tile);

            ground = compute_diagnostic@READ_STATVAR_FROM_OUT(ground, tile);

        end
        
        function ground = check_trigger(ground, tile)
            ground = check_trigger@READ_STATVAR_FROM_OUT(ground, tile);
            ground.BGC = check_trigger(ground.BGC, tile);
        end
        

    end
    
end
