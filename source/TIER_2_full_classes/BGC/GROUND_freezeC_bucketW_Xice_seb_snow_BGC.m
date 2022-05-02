%========================================================================
% CryoGrid GROUND class GROUND_freezeC_bucketW_Xice_seb_snow_BGC
% enables coupling to biogeochemistry BGC classes
% S. Westermann, November 2021
%========================================================================

classdef GROUND_freezeC_bucketW_Xice_seb_snow_BGC < GROUND_freezeC_bucketW_Xice_seb_snow

    properties
        BGC
        IA_BGC
    end
    
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

        
        function ground = provide_PARA(ground)
            
            ground = provide_PARA@GROUND_freezeC_bucketW_Xice_seb_snow(ground);
            
            ground.PARA.BGC_CLASS = [];
            ground.PARA.target_grid_cell_size = 0.05;
        end
        
        function ground = provide_STATVAR(ground)
            
            ground = provide_STATVAR@GROUND_freezeC_bucketW_Xice_seb_snow(ground);

        end
        
        function ground = provide_CONST(ground)
            
            ground = provide_CONST@GROUND_freezeC_bucketW_Xice_seb_snow(ground);
            
        end
        
        function ground = finalize_init(ground, tile) 

            ground = finalize_init@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);

            class_handle = str2func(ground.PARA.BGC_CLASS);
            ground.BGC = class_handle(); 
            ground.BGC.PARENT = ground;
            
            ground.BGC = provide_PARA(ground.BGC);
            ground.BGC = provide_STATVAR(ground.BGC);
            ground.BGC = provide_CONST(ground.BGC);
            ground.BGC = finalize_init(ground.BGC, tile);
            
            ground.IA_BGC = IA_BGC_Xice();
            ground.IA_BGC.BGC = ground.BGC;
            ground.IA_BGC.GROUND = ground;
            ground.BGC.IA_BGC = ground.IA_BGC;
            finalize_init(ground.IA_BGC, tile);
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            ground = get_boundary_condition_u@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            ground.BGC = get_boundary_condition_u(ground.BGC, tile);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            ground = get_boundary_condition_l@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
           ground.BGC = get_boundary_condition_l(ground.BGC, tile); 
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground = get_derivatives_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
           ground.BGC = get_derivatives_prognostic(ground.BGC, tile); 
        end
        
        function timestep = get_timestep(ground, tile) 
            timestep_ground = get_timestep@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            timestep_BGC = get_timestep(ground.BGC, tile);
            timestep = min(timestep_ground, timestep_BGC);
        end
        
        function ground = advance_prognostic(ground, tile)    
            ground = advance_prognostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            ground.BGC = advance_prognostic(ground.BGC, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            ground = compute_diagnostic_first_cell@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            ground.BGC = compute_diagnostic_first_cell(ground.BGC, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)
            
            ground.BGC = compute_diagnostic(ground.BGC, tile);
                        
            ground = compute_diagnostic@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);

        end
        
        function ground = check_trigger(ground, tile)
            ground = check_trigger@GROUND_freezeC_bucketW_Xice_seb_snow(ground, tile);
            ground.BGC = check_trigger(ground.BGC, tile);
        end
    
        function ground = reset_timestamps(ground, tile) %used e.g. with TILE_BUILDER update_forcing_out
            ground.BGC = reset_time(ground.BGC, tile);
        end
        
        
    end
    
end
