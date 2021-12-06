%========================================================================
% CryoGrid GROUND class GROUND_store_flip_flop_singleClass_BGC
% reads variables from FLIP_FLOP files, to be used in conjuction with
% ..._FLIP_FLOP classes, from where it is called. It is used with
% biogeochemeistry classes.
% S. Westermann, November 2021
%========================================================================

classdef GROUND_store_flip_flop_singleClass_BGC < GROUND_store_flip_flop_singleClass

    properties
        BGC
        IA_BGC
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
         function ground = provide_PARA(ground)
            
            ground = provide_PARA@GROUND_store_flip_flop_singleClass(ground);

            ground.PARA.BGC_CLASS = [];
            
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@GROUND_store_flip_flop_singleClass(ground);
        end
        
        function ground = provide_CONST(ground)
            
            ground.CONST.day_sec = [];
        end
        
        function ground = finalize_init(ground, tile) 
             ground = finalize_init@GROUND_store_flip_flop_singleClass(ground, tile);

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            ground = get_boundary_condition_u@GROUND_store_flip_flop_singleClass(ground, tile);
            ground.BGC = get_boundary_condition_u(ground.BGC, tile);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_store_flip_flop_singleClass(ground, tile);
            ground.BGC = get_boundary_condition_l(ground.BGC, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivatives_prognostic@GROUND_store_flip_flop_singleClass(ground, tile);
            ground.BGC = get_derivatives_prognostic(ground.BGC, tile);
        end
        
        function timestep = get_timestep(ground, tile)
            timestep_ground = get_timestep@GROUND_store_flip_flop_singleClass(ground, tile);
            timestep_BGC = get_timestep(ground.BGC, tile);
            timestep = min(timestep_ground, timestep_BGC);
        end
        
        function ground = advance_prognostic(ground, tile) 
            ground = advance_prognostic@GROUND_store_flip_flop_singleClass(ground, tile);
            ground.BGC = advance_prognostic(ground.BGC, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            ground = compute_diagnostic_first_cell@GROUND_store_flip_flop_singleClass(ground, tile);
            ground.BGC = compute_diagnostic_first_cell(ground.BGC, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)

            ground.BGC = compute_diagnostic(ground.BGC, tile);
            ground = compute_diagnostic@GROUND_store_flip_flop_singleClass(ground, tile);

        end
        
        function ground = check_trigger(ground, tile)

            ground.BGC = check_trigger(ground.BGC, tile);
            ground = switch2model_flip_flop_BGC(ground, tile);
        end

    end
    
end
