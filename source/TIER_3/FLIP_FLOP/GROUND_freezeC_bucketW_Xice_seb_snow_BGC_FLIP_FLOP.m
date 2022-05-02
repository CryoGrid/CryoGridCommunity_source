%========================================================================
% CryoGrid TIER 3 GROUND class combining snow and biogeochemistry coupling
% with write/read capabilities of the FLIP_FLOP Tier 1 class
% This class can serve as a template to enable FLIP_FLOP also for other
% GROUND classes
% S. Westermann, November 2021
%========================================================================

classdef GROUND_freezeC_bucketW_Xice_seb_snow_BGC_FLIP_FLOP < GROUND_freezeC_bucketW_Xice_seb_snow_BGC & FLIP_FLOP

    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            ground = provide_PARA@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground);
            ground = provide_PARA_flip_flop(ground);
        end
        
        function ground = provide_CONST(ground)
            ground = provide_CONST@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground);
            ground = provide_CONST_flip_flop(ground);
        end
        
        function ground = provide_STATVAR(ground)
            ground = provide_STATVAR@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground);
        end
        
        function ground = finalize_init(ground, tile)
            ground = finalize_init@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
            ground = finalize_init_flip_flop(ground, tile);

        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
             ground = get_boundary_condition_u@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            [ground, S_up] = penetrate_SW@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, S_down);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)

            ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile); %call normal function

        end
        
        function timestep = get_timestep(ground, tile)

            timestep = get_timestep@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
            timestep = min(timestep, get_timestep_flip_flop(ground, tile));

        end
        
        function ground = advance_prognostic(ground, tile)

            ground =  advance_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);

        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            
            ground = compute_diagnostic_first_cell@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
        
        end
        
        function ground = compute_diagnostic(ground, tile)

            ground = compute_diagnostic@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);

        end
        
        function ground = check_trigger(ground, tile)
            
            ground = check_trigger@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
            
            ground = store_flip_flop(ground, tile);
            ground = write_store_flip_flop(ground, tile);
            ground = switch2read_STATVAR_flip_flop_BGC(ground, tile);
        end
        
        function ground = reset_timestamps(ground, tile) %used e.g. with TILE_BUILDER update_forcing_out
            ground = reset_timestamps@GROUND_freezeC_bucketW_Xice_seb_snow_BGC(ground, tile);
            ground = finalize_init_flip_flop(ground, tile);            
        end
        
    end
end