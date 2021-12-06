%========================================================================
% CryoGrid GROUND class 

% S. Westermann, October 2020
%========================================================================

classdef GROUND_freezeC_bucketW_Xice_seb_snow_FLIP_FLOP < GROUND_freezeC_bucketW_Xice_seb_snow & FLIP_FLOP

    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            ground = provide_PARA@GROUND_freezeC_bucketW_Xice_seb_snow(ground);
            ground = provide_PARA_flip_flop(ground);
        end
        
        function ground = provide_CONST(ground)
            ground = provide_CONST@GROUND_freezeC_bucketW_Xice_seb_snow(ground);
            ground = provide_CONST_flip_flop(ground);
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@GROUND_freezeC_bucketW_Xice_seb_snow(ground);
        end
        
        function ground = finalize_init(ground, tile)
            ground = finalize_init@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            ground = finalize_init_flip_flop(ground, tile);

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
             ground = get_boundary_condition_u@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@GROUND_freezeC_bucketW_Xice_seb_snow(ground, S_down);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)

            ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile); %call normal function

        end
        
        function timestep = get_timestep(ground, tile)

            timestep = get_timestep@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            timestep = min(timestep, get_timestep_flip_flop(ground, tile));

        end
        
        function ground = advance_prognostic(ground, tile)

            ground =  advance_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);

        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            
            ground = compute_diagnostic_first_cell@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
        
        end
        
        function ground = compute_diagnostic(ground, tile)

            ground = compute_diagnostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);

        end
        
        function ground = check_trigger(ground, tile)
            
            ground = check_trigger@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            
            ground = store_flip_flop(ground, tile);
            ground = write_store_flip_flop(ground, tile);
            ground = switch2read_STATVAR_flip_flop(ground, tile);
        end
        
    end
end