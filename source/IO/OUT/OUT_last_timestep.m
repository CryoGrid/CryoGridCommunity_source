%========================================================================
% CryoGrid OUT class OUT_last_timestep
% can be used to periodically store the entire CryoGrid stratigraphy, with 
% the ouput file getting overwritten each time. This is useful for being 
% able to recover and restart simulations (e.g. when the cluster is shut 
% down due to maintenances) in case of a long runtime. It is also posisble 
% to only store the final state of the CryoGrid stratigraphy, e.g. to run
% ensembles starting from the final state of an existing TILE class. Use
% TILE_BUILDER class "restart_OUT_last_timestep" to restart simulations.
% S. Westermann, October 2020
%========================================================================


classdef OUT_last_timestep < matlab.mixin.Copyable

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
%         function out = initialize_excel(out)
%             
%         end
        
        
        function out = provide_PARA(out)         

            out.PARA.save_timestep = []; %if empty save final state at the end of the run, so that it can serve as initial condition for new runs
            out.PARA.tag = [];

        end
		
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)
		
			forcing = tile.FORCING;
            if isempty(out.PARA.save_timestep) || isnan(out.PARA.save_timestep)
                out.OUTPUT_TIME = forcing.PARA.end_time;
            else
                out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.save_timestep;
            end
            out.SAVE_TIME = forcing.PARA.end_time;

        end
        
        %-------time integration----------------
        
        function out = store_OUT(out, tile)
            t = tile.t;
            out_tag = out.PARA.tag;
            
            if t==out.SAVE_TIME || t == out.OUTPUT_TIME
                
                disp([datestr(t)])
                
                run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
                result_path = tile.PARA.result_path;
                
                
                out.STRATIGRAPHY = copy(tile);


                if ~(exist([result_path run_name])==7)
                    mkdir([result_path result_path])
                end
                if isempty(out_tag) || all(isnan(out_tag))
                    save([result_path run_name '/' run_name '_last_timestep.mat'], 'out')
                else
                    save([result_path run_name '/' run_name '_' out_tag '_last_timestep.mat'], 'out')
                end
                
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.save_timestep;
            end
            
        end
        
                %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
           
            out.PARA.default_value.save_timestep = {''};
            out.PARA.comment.save_timestep = {'in days, if empty save final state at the end of the run, so that it can serve as initial condition for new runs'};
            
            out.PARA.default_value.tag = {''};
            out.PARA.comment.tag = {'additional tag added to file name'};
        end

        
    end
end