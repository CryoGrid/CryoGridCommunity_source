%========================================================================
% CryoGrid OUT class
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

        
    end
end