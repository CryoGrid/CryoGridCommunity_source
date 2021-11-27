%========================================================================
% CryoGrid GROUND class GROUND_store_flip_flop_singleClass
% reads variables from FLIP_FLOP files, to be used in conjuction with
% ..._FLIP_FLOP classes, from where it is called
% S. Westermann, November 2021
%========================================================================

classdef GROUND_store_flip_flop_singleClass < FLIP_FLOP


    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------

         function ground = provide_PARA(ground)
            

        end
        
        function ground = provide_STATVAR(ground)

        end
        
        function ground = provide_CONST(ground)

        end
        
        function ground = finalize_init(ground, tile) 
        
              ground.TEMP.next_read_time = tile.t + ground.PARA.store_interval;
              ground.TEMP.next_file_time = datenum([ground.PARA.save_date num2str(str2num(datestr(tile.t, 'yyyy')) + 1)], 'dd.mm.yyyy');
              ground.TEMP.year_count = 1;
              run_name = tile.PARA.run_name;
              result_path = tile.PARA.result_path;
              load([result_path run_name '/' run_name '_STORE_' num2str(ground.TEMP.year_count) '.mat']);
              ground.STORE = store;
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            
        end
        
%         function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
%             S_up = 0;
%             %[ground, S_up] = penetrate_SW_no_transmission(ground, S_down);
%         end
        
        function ground = get_boundary_condition_l(ground, tile)

        end
        
        function ground = get_derivatives_prognostic(ground, tile)
            
        end
        
        function timestep = get_timestep(ground, tile)

            timestep = get_timestep_flip_flop2(ground, tile);

        end
        
        function ground = advance_prognostic(ground, tile)
            
        end
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            
        end
        
        function ground = compute_diagnostic(ground, tile)
            
            if tile.t == ground.TEMP.next_read_time %read and replace STATVAR
                offset_years = floor((tile.t - ground.STORE.stratigraphy{1,end}.STATVAR.timestamp) ./ 365.25)+1;
                index = round((tile.t - offset_years .* 365.25 - ground.STORE.stratigraphy{1,1}.STATVAR.timestamp) ./ ground.PARA.store_interval) + 1;
                index(index<=0) = size(ground.STORE.stratigraphy,2) - index;
                index(index>size(ground.STORE.stratigraphy,2)) = index(index>size(ground.STORE.stratigraphy,2)) - floor(index./size(ground.STORE.stratigraphy,2)) .* size(ground.STORE.stratigraphy,2);
             
                ground.STATVAR = ground.STORE.stratigraphy{1,index}.STATVAR; 
                
                ground.TEMP.next_read_time = tile.t + ground.PARA.store_interval;
                
            end
            if tile.t>=ground.TEMP.next_file_time %read next file
                ground.TEMP.year_count = ground.TEMP.year_count + 1;
                if ground.TEMP.year_count> ground.PARA.number_of_model_years - ground.PARA.start_store_year
                    ground.TEMP.year_count = 1;
                end
                run_name = tile.PARA.run_name;
                result_path = tile.PARA.result_path;
                load([result_path run_name '/' run_name '_STORE_' num2str(ground.TEMP.year_count) '.mat']);
                ground.STORE = store;
                
                save_day = str2num(ground.PARA.save_date(1:2));
                save_month = str2num(ground.PARA.save_date(4:5));
                [year, ~,~] = datevec(ground.TEMP.next_file_time);
                ground.TEMP.next_file_time = datenum(year+1, save_day, save_month);
               % ground.TEMP.next_file_time = datenum([ground.PARA.save_date num2str(str2num(datestr(ground.TEMP.next_file_time, 'yyyy')) + 1)], 'dd.mm.yyyy');
                
            end
                
        end
        
        function ground = check_trigger(ground, tile)
            ground = switch2model_flip_flop(ground, tile);
        end
        

    end
    
end
