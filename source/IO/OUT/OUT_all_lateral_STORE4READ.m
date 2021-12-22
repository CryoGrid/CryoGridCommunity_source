%========================================================================
% CryoGrid OUT class defining storage format of the output 
% OUT_all_lateral_STORE4READ stores identical copies of all GROUND classses 
% (including STATVAR, TEMP, PARA) in the stratigraphy for each output 
% timestep and copies of all LATERAL classes.
% Other than that, it is identical to OUT_all.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, June 2021
%========================================================================


classdef OUT_all_lateral_STORE4READ < matlab.mixin.Copyable

    properties
		out_index
        STRATIGRAPHY
        LATERAL
        TIMESTAMP
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
		
	end
    
    
    methods
		

        %-------initialization--------------
        
% 		function out = initialize(out)
%             % INITIALIZE  Initializes all properties needed by the class.
% 			
% 			out.STRATIGRAPHY = [];
%             out.LATERAL=[];
% 			out.TIMESTAMP = [];
% 			out.OUTPUT_TIME = [];
% 			out.SAVE_TIME = [];
%             out = out.initialize_PARA();
%             out = out.initialize_TEMP();
%         end
        
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.number_of_skipped_classes = 1;
            out.PARA.tag = [];
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end

		function  out = finalize_init(out, tile)

            forcing = tile.FORCING;

			out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = tile.FORCING.PARA.end_time;
            else
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.FORCING.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            out.TEMP = struct();
        end
        
        %-------time integration----------------
		
		function out = store_OUT(out, tile)
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             forcing = tile.FORCING;
             run_name = tile.PARA.run_name;
             result_path = tile.PARA.result_path;
             timestep = tile.timestep;
             out_tag = out.PARA.tag;

            
            if t==out.OUTPUT_TIME
                % It is time to collect output
                % Store the current state of the model in the out structure.

                disp([datestr(t)])
                out.TIMESTAMP=[out.TIMESTAMP t];
                
                new_BOTTOM = Bottom();
                CURRENT=BOTTOM.PREVIOUS;
                new_CURRENT = new_BOTTOM;
                count = 0;
                while ~isequal(CURRENT, TOP)
                    if count >= out.PARA.number_of_skipped_classes
                        new_CURRENT.PREVIOUS = copy(CURRENT);
                        new_CURRENT.PREVIOUS.NEXT = new_CURRENT;
                        new_CURRENT = new_CURRENT.PREVIOUS;
                        if isprop(new_CURRENT, 'LUT')
                            new_CURRENT.LUT =[];  %remove look-up tables, runs out of memory otherwise
                        end
                        new_CURRENT.PREVIOUS = [];
                        new_CURRENT.IA_PREVIOUS =[];
                        new_CURRENT.IA_NEXT = [];
                    end
                    count = count + 1;
                    CURRENT = CURRENT.PREVIOUS;
                end
                new_TOP = copy(TOP);
                new_TOP.STORE = [];
                new_TOP.LATERAL = [];
                new_TOP.FORCING = [];
                new_CURRENT.PREVIOUS = new_TOP;
                new_TOP.NEXT = new_CURRENT;
                
                out.STRATIGRAPHY{1,size(out.STRATIGRAPHY,2)+1} = new_BOTTOM;
                
                %lateral, read only STATVAR and PARA---
                result={};
                ia_classes=TOP.LATERAL.IA_CLASSES;
                for i=1:size(ia_classes,1)
                    res = copy(ia_classes{i,1});
                    vars = fieldnames(res);
                    for j=1:size(vars,1)
                        if ~strcmp(vars{j,1}, 'PARA') && ~strcmp(vars{j,1}, 'STATVAR')
                           res.(vars{j,1}) = []; 
                        end
                    end
                    result=[result; {res}];                    
                end
                out.LATERAL{1,size(out.LATERAL, 2)+1} = result;
                %---
                
				% Set the next OUTPUT_TIME
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
				
                if t>=out.SAVE_TIME 
					% It is time to save all the collected model output to disk
                    if ~(exist([result_path run_name])==7)
                        mkdir([result_path run_name])
                    end
                    if isempty(out_tag) || all(isnan(out_tag))
                        save([result_path run_name '/' run_name '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    else
                        save([result_path run_name '/' run_name '_' out_tag '_' datestr(t,'yyyymmdd') '.mat'], 'out')
                    end
                    
					
					% Clear the out structure
                    out.STRATIGRAPHY=[];
                    out.LATERAL=[];
                    out.TIMESTAMP=[];
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