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

            out.PARA.save_timestep = [];

        end
		
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)
		
			forcing = tile.FORCING;

			out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.save_timestep;
            out.SAVE_TIME = forcing.PARA.end_time;

        end
        
        %-------time integration----------------
        
        function out = store_OUT(out, tile)
            t = tile.t;
            if t==out.SAVE_TIME || t == out.OUTPUT_TIME
                
                disp([datestr(t)])
                
                run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
                result_path = tile.PARA.result_path;
                
                
                out.STRATIGRAPHY = copy(tile);


                if ~(exist([result_path run_name])==7)
                    mkdir([result_path result_path])
                end
                save([result_path run_name '/' run_name '_last_timestep.mat'], 'out')
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.save_timestep;
            end
            
        end

        
    end
end