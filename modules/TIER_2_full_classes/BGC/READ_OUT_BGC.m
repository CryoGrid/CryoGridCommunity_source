%========================================================================
% CryoGrid GROUND class READ_OUT_BGC
% heat conduction, bucket water scheme, freeze curve based on
% freezing=drying assumption, surface energy balance, excess ice
% testing coupling to biogeochemistry (BGC) classes
% S. Westermann, October 2020
%========================================================================

classdef READ_OUT_BGC < INITIALIZE

    properties
        READ_OUT
        RUN_PARA
        RUN_CONST
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        function ground = READ_OUT_BGC(index, pprovider, cprovider, forcing)  
            ground@INITIALIZE(index, pprovider, cprovider, forcing);
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
            %ground.PARA.out_number_of_skipped_classes = [];
            
            %only needs compute_diagnostic: read the stratigraphy from OUT and replace the classes in the
            %present stratigraphy bottom-up; other than that, call all the functions
            %for the BGC class
        end
        
        function ground = provide_STATVAR(ground)
            ground.STATVAR.dummy = [];
        end
        
        function ground = provide_CONST(ground)
            ground.CONST.day_sec = [];
        end
        
        function ground = finalize_init(ground, forcing) 
            ground.STATVAR.year_list = [ground.PARA.start_year:ground.PARA.out_save_interval:ground.PARA.end_year];
            ground.STATVAR.year_index = 1; 
            
            filename = [ground.PARA.run_number '_' datestr(datenum([ground.PARA.out_save_date num2str(ground.STATVAR.year_list(ground.STATVAR.year_index))], 'dd.mm.yyyy'), 'yyyymmdd') '.mat'];
            load([ground.PARA.result_path '/' ground.PARA.run_number '/' filename]);
            ground.READ_OUT = out;  %load the first file
            
            ground.PARA.time_offset = datenum([ground.PARA.out_save_date '2002'], 'dd.mm.yyyy')-datenum(2000,1,1);
            
            ground.RUN_PARA = ground.PARA;
            ground.RUN_CONST = ground.CONST;
            ground.RUN_PARA.next_read_time = forcing.PARA.start_time + ground.RUN_PARA.timestep; 
            ground.RUN_PARA.active = 1;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, forcing)

        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            S_up = 0;
            %[ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, forcing)

        end
        
        function ground = get_derivatives_prognostic(ground)

        end
        
        function timestep = get_timestep(ground) 
            if ground.RUN_PARA.active
                timestep = ground.RUN_PARA.next_read_time - ground.PREVIOUS.TIME;  %Refers to the field TOP.TIME
                timestep = timestep .* ground.RUN_CONST.day_sec;
            else
               timestep = Inf;
           end
        end
        
        function ground = advance_prognostic(ground, timestep) 

        end
        
        function ground = compute_diagnostic_first_cell(ground, forcing)

        end
       
        function ground = compute_diagnostic(ground, forcing)
            %[year, month, day, hour, minute] = datevec(ground.PREVIOUS.TIME - ground.PARA.time_offset);
            time_adjusted = (ground.PREVIOUS.TIME - ground.RUN_PARA.time_offset); %adjust this!!
            [year, month, day] = datevec(time_adjusted);
            doy_since_end =datenum(year+1, 1, 1) - time_adjusted;
            
            i= size(ground.READ_OUT.TIMESTAMP,2) -  doy_since_end./ground.RUN_PARA.out_output_timestep

            ground.RUN_PARA.next_read_time = ground.RUN_PARA.next_read_time + ground.RUN_PARA.timestep; 
            
            CURRENT = ground;
            while ~isequal(class(CURRENT), 'Top')
                CURRENT = CURRENT.PREVIOUS;
            end
            old_Top = CURRENT;
            
            CURRENT = ground.READ_OUT.STRATIGRAPHY{1,i}.PREVIOUS;
            
            ground.STATVAR = CURRENT.STATVAR;
            ground.PARA = CURRENT.PARA;
            ground.TEMP = CURRENT.TEMP;
            ground.CONST = CURRENT.CONST;
            
            while ~isequal(class(CURRENT.PREVIOUS), 'Top')
                CURRENT = CURRENT.PREVIOUS;
                new_GROUND = READ_OUT_BGC(-1,0,0,0);
                new_GROUND.STATVAR = CURRENT.STATVAR;
                new_GROUND.PARA = CURRENT.PARA;
                new_GROUND.TEMP = CURRENT.TEMP;
                new_GROUND.CONST = CURRENT.CONST;
                new_GROUND.RUN_PARA.active = 0;
                new_GROUND.NEXT = ground;
                new_GROUND.IA_NEXT = IA_DO_NOTHING();
                ground.PREVIOUS = new_GROUND;
                ground.IA_PREVIOUS = IA_DO_NOTHING();
                ground = new_GROUND;
            end
            ground.PREVIOUS = old_Top;
            old_Top.NEXT = ground;
            

        end
        
        function ground = check_trigger(ground, forcing)

        end
        

    end
    
end
