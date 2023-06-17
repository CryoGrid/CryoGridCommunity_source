%========================================================================
% CryoGrid OUT class OUT_TwaterIce defining model output and storage 
% OUT_TwaterIce interpolates the desired output variables to a fixed 
% vertical grid. The output files are in Matlab (".mat") format.
%
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, June 2021
% S. Westermann, Oct 2022
%========================================================================


classdef OUT_SEB < matlab.mixin.Copyable
 

    properties
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        result
        timestamp
        
    end
    
    
    methods
                
        function out = provide_PARA(out)         

            out.PARA.variables = [];
            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.tag = [];
        end

        
        function out = provide_CONST(out)

        end

        
        function out = provide_STATVAR(out)

        end

        
        function out = finalize_init(out, tile)

            forcing = tile.FORCING;
            
            % Set the next (first) output time. This is the next (first) time output
            % is collected (in memory) for later storage to disk.
            out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            %make struct "result" and initialize all variables defined by
            %the user as empty arrays
            out.timestamp = [];
            for i =  1:size(out.PARA.variables,1)
               out.result.(out.PARA.variables{i,1}) = []; 
            end
            
            % Set the next (first) save time. This is the next (first) time all the
            % collected output is saved to disk.
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
                        
        end
        
        %---------------time integration-------------
        
%         function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
            t = tile.t;
            TOP = tile.TOP; 
            BOTTOM = tile.BOTTOM;
            forcing = tile.FORCING;
            run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
            result_path = tile.PARA.result_path;            
            timestep = tile.timestep;
            out_tag = out.PARA.tag;
            
            if t>=out.OUTPUT_TIME
                % It is time to collect output
                % Store the current state of the model in the out structure.

%                 disp([datestr(t)])                
                out.timestamp = [out.timestamp t];
                
                CURRENT =TOP.NEXT;
                for i=1:size(out.PARA.variables,1)
                    if strcmp(out.PARA.variables{i,1}, 'Sin') || strcmp(out.PARA.variables{i,1}, 'Lin')
                       out.result.(out.PARA.variables{i,1}) = [out.result.(out.PARA.variables{i,1}) tile.FORCING.TEMP.(out.PARA.variables{i,1})];
                    else
                        out.result.(out.PARA.variables{i,1}) = [out.result.(out.PARA.variables{i,1}) CURRENT.STATVAR.(out.PARA.variables{i,1})];
                    end
                end
                
                % Set the next OUTPUT_TIME
                out.OUTPUT_TIME = min(out.SAVE_TIME, out.OUTPUT_TIME + out.PARA.output_timestep);
                
                if t>=out.SAVE_TIME
                    % It is time to save all the collected model output to disk
                     
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    OUT_SEB = out.result;
                    OUT_SEB.timestamp = out.timestamp;
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_SEB_' datestr(t,'yyyymmdd') '.mat'], 'OUT_SEB')
                    else
                        save([result_path run_name '/' run_name '_SEB_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'OUT_SEB')
                    end
                    
                    % Clear the out structure
                    out.MISC=[];
                    %make struct "result" and initialize all variables defined by
                    %the user as empty arrays
                    out.timestamp = [];
                    for i =  1:size(out.PARA.variables,1)
                        out.result.(out.PARA.variables{i,1}) = [];
                    end
                    out.result.depths = [];
                    out.result.class_number = [];
                    
                    if ~isnan(out.PARA.save_interval)
                        % If save_interval is defined, uptate SAVE_TIME for next save opertion 
                        out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(out.SAVE_TIME,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                        % If save_interval is not defined, we will save at the very end of the model run
                        % and thus do not need to update SAVE_TIME (update would fail because save_interval is nan)
					end
                end
            end
        end
        
                

        %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
            
            out.PARA.comment.variables = {'select output variables: Sin, Sout, Lin, Lout, Qe, Qh are supported'};
            out.PARA.options.variables.name = 'H_LIST';
            out.PARA.options.variables.entries_x = {'Sin' 'Sout' 'Lin' 'Lout' 'Qh' 'Qe'};
           
            out.PARA.default_value.output_timestep = {0.25};
            out.PARA.comment.output_timestep = {'timestep of output [days]'};

            out.PARA.default_value.save_date = {'01.09.'};
            out.PARA.comment.save_date = {'date (dd.mm.) when output file is written'};
            
            out.PARA.default_value.save_interval = {1};
            out.PARA.comment.save_interval = {'interval of output files [years]'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end
        

    end
end