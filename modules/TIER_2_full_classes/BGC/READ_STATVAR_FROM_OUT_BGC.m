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
        
        function ground = READ_STATVAR_FROM_OUT_BGC(index, pprovider, cprovider, forcing)  
            ground@READ_STATVAR_FROM_OUT(index, pprovider, cprovider, forcing);
        end
        
         function ground = provide_PARA(ground)
            
            ground.PARA.start_year = [];
            ground.PARA.end_year = [];
            ground.PARA.timestep = []; %must  be multiple of output timestep
            ground.PARA.run_number = []; %run _number to read - if empty, use own run number (only works if this run 
            ground.PARA.result_path = [];
            
            ground.PARA.out_output_timestep = [];
            ground.PARA.out_save_date = [];
            ground.PARA.out_save_interval = [];
            
            ground.PARA.BGC_CLASS = [];
            
        end
        
        function ground = provide_STATVAR(ground)
            ground.STATVAR.dummy = [];
        end
        
        function ground = provide_CONST(ground)
            ground.CONST.day_sec = [];
        end
        
        function ground = finalize_init(ground, forcing) 
            ground = finalize_init@READ_STATVAR_FROM_OUT(ground, forcing);

            class_handle = str2func(ground.RUN_PARA.BGC_CLASS);
            ground.BGC = class_handle(-1,0,0,0); 
            ground.BGC.PARENT = ground;
            %remove this in the end
            ground.BGC = provide_PARA(ground.BGC);
            ground.BGC = provide_STATVAR(ground.BGC);
            ground.BGC = provide_CONST(ground.BGC);
            ground.BGC = finalize_init(ground.BGC, forcing);
            
            ground.IA_BGC = IA_BGC_simple();
            ground.IA_BGC.BGC = ground.BGC;
            ground.IA_BGC.GROUND = ground;
            ground.BGC.IA_BGC = ground.IA_BGC;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
           
            ground.BGC = get_boundary_condition_u(ground.BGC, tile);
        end
        
        
        function ground = get_boundary_condition_l(ground, tile)
            
            ground.BGC = get_boundary_condition_l(ground.BGC, tile);
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            ground.BGC = get_derivatives_prognostic(ground.BGC, tile);
        end
        
        function timestep = get_timestep(ground, tile)
            timestep_phyiscal = get_timestep@READ_STATVAR_FROM_OUT(ground, tile);
            timestep_BGC = get_timestep(ground.BGC, tile);
            timestep = min(timestep_phyiscal, timestep_BGC);
        end
        
        function ground = advance_prognostic(ground, tile) 
            ground.BGC = advance_prognostic(ground.BGC, tile);
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            ground.BGC = compute_diagnostic_first_cell(ground.BGC, tile);
        end
       
        function ground = compute_diagnostic(ground, tile)

            ground.BGC = compute_diagnostic(ground.BGC, tile);

                        
            ground = compute_diagnostic@READ_STATVAR_FROM_OUT(ground, tile);


        end
        
        function ground = check_trigger(ground, tile)
            ground.BGC = check_trigger(ground.BGC, tile);
        end
        

    end
    
end
