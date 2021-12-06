%========================================================================
% CryoGrid GROUND class READ_STATVAR_FROM_OUT
% reads variables from OUT files of the OUT class OUT_all_lateral_STORE4READ
% S. Westermann, November 2020
%========================================================================

classdef READ_STATVAR_FROM_OUT < BASE

    properties
        RUN_PARA
        RUN_CONST
        READ_OUT
    end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        

        
         function ground = provide_PARA(ground)
            
            ground.PARA.start_year = [];
            ground.PARA.end_year = [];
            
            ground.PARA.timestep = []; %must  be multiple of output timestep
            ground.PARA.run_name = []; %run _number to read - if empty, use own run number (only works if this run 
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
        
        function ground = finalize_init(ground, tile) 
        
            if isempty(ground.PARA.result_path) || isnan(ground.PARA.result_path(1,1))
                ground.PARA.result_path = tile.PARA.result_path;
            end
            if isempty(ground.PARA.run_name) || isnan(ground.PARA.run_name(1,1))
                ground.PARA.run_name = tile.PARA.run_name;
            end
            
            ground.PARA.year_list = [ground.PARA.start_year:ground.PARA.out_save_interval:ground.PARA.end_year];
            ground.PARA.year_index = 1; 

            filename = [ground.PARA.run_name '_' datestr(datenum([ground.PARA.out_save_date num2str(ground.PARA.year_list(ground.PARA.year_index))], 'dd.mm.yyyy'), 'yyyymmdd') '.mat'];
            load([ground.PARA.result_path ground.PARA.run_name '/' filename]);
            ground.READ_OUT = out;  %load the first file

            ground.PARA.time_offset = datenum([ground.PARA.out_save_date num2str(ground.PARA.year_list(ground.PARA.year_index))], 'dd.mm.yyyy') - datenum(ground.PARA.year_list(ground.PARA.year_index),1,1);
            ground.RUN_PARA = ground.PARA;
            ground.RUN_CONST = ground.CONST;
            ground.RUN_PARA.next_read_time = tile.FORCING.PARA.start_time + ground.RUN_PARA.timestep; 
            ground.RUN_PARA.active = 1;
            
            %assign the first stratigraphy
            CURRENT = ground.READ_OUT.STRATIGRAPHY{1,1}.PREVIOUS;
            ground.STATVAR = CURRENT.STATVAR;
            ground.PARA = CURRENT.PARA;
            ground.TEMP = CURRENT.TEMP;
            ground.CONST = CURRENT.CONST;
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            
        end
        
        function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
            S_up = 0;
            %[ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
        end
        
        function ground = get_boundary_condition_l(ground, tile)
            
        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            
        end
        
        function timestep = get_timestep(ground, tile)
            if ground.RUN_PARA.active
                timestep = ground.RUN_PARA.next_read_time - tile.t;  %Refers to the field TOP.TIME
                timestep = timestep .* ground.RUN_CONST.day_sec;
            else
                timestep = Inf;
            end
        end
        
        function ground = advance_prognostic(ground, tile)
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            
        end
        
        function ground = compute_diagnostic(ground, tile)
            if ground.RUN_PARA.active
                index = floor((tile.t - tile.FORCING.PARA.start_time) ./ (365.25 .* ground.RUN_PARA.out_save_interval));
                index = mod(index, size(ground.RUN_PARA.year_list,2)) + 1;
                
                %load new out file
                if ground.RUN_PARA.year_index ~= index
                    ground.RUN_PARA.year_index = index;
                    filename = [ground.RUN_PARA.run_name '_' datestr(datenum([ground.RUN_PARA.out_save_date num2str(ground.RUN_PARA.year_list(ground.RUN_PARA.year_index))], 'dd.mm.yyyy'), 'yyyymmdd') '.mat'];
                    disp(['reading ' filename])
                    load([ground.RUN_PARA.result_path '/' ground.RUN_PARA.run_name '/' filename]);
                    ground.READ_OUT = out;
                    ground.RUN_PARA.time_offset = datenum([ground.RUN_PARA.out_save_date num2str(ground.RUN_PARA.year_list(ground.RUN_PARA.year_index))], 'dd.mm.yyyy') - datenum(ground.RUN_PARA.year_list(ground.RUN_PARA.year_index),1,1);
                end
                
                time_adjusted = (tile.t - ground.RUN_PARA.time_offset);
                [year, month, day] = datevec(time_adjusted);
                doy_since_end =datenum(year+1, 1, 1) - time_adjusted;
                
                i= mod(round(size(ground.READ_OUT.TIMESTAMP,2) -  doy_since_end./ground.RUN_PARA.out_output_timestep) - 1, size(ground.READ_OUT.TIMESTAMP,2)) + 1;
                
                
                ground.RUN_PARA.next_read_time = ground.RUN_PARA.next_read_time + ground.RUN_PARA.timestep;
                
                %find TOP class
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
                
                %new from here
                CURRENT2 = ground;
                
                while ~isequal(class(CURRENT.PREVIOUS), 'Top')
                    CURRENT = CURRENT.PREVIOUS;
                    new_GROUND = READ_STATVAR_FROM_OUT();
                    new_GROUND.STATVAR = CURRENT.STATVAR;
                    new_GROUND.PARA = CURRENT.PARA;
                    %new_GROUND.RUN_PARA = ground.RUN_PARA;
                    new_GROUND.TEMP = CURRENT.TEMP;
                    new_GROUND.CONST = CURRENT.CONST;
                    new_GROUND.RUN_PARA.active = 0;
                    
                    new_GROUND.NEXT = CURRENT2;
                    new_GROUND.IA_NEXT = IA_DO_NOTHING();
                    CURRENT2.PREVIOUS = new_GROUND;
                    CURRENT2.IA_PREVIOUS = IA_DO_NOTHING();
                    CURRENT2 = new_GROUND;
                end
                CURRENT2.PREVIOUS = old_Top;
                old_Top.NEXT = CURRENT2;
                
%                 while ~isequal(class(CURRENT.PREVIOUS), 'Top')
%                     CURRENT = CURRENT.PREVIOUS;
%                     new_GROUND = READ_STATVAR_FROM_OUT();
%                     new_GROUND.STATVAR = CURRENT.STATVAR;
%                     new_GROUND.PARA = CURRENT.PARA;
%                     %new_GROUND.RUN_PARA = ground.RUN_PARA;
%                     new_GROUND.TEMP = CURRENT.TEMP;
%                     new_GROUND.CONST = CURRENT.CONST;
%                     new_GROUND.RUN_PARA.active = 0;
%                     new_GROUND.NEXT = ground;
%                     new_GROUND.IA_NEXT = IA_DO_NOTHING();
%                     ground.PREVIOUS = new_GROUND;
%                     ground.IA_PREVIOUS = IA_DO_NOTHING();
%                     ground = new_GROUND;
%                 end
%                 ground.PREVIOUS = old_Top;
%                 old_Top.NEXT = ground;
                
            end

        end
        
        function ground = check_trigger(ground, tile)

        end
        

    end
    
end
