%========================================================================
% CryoGrid OUT class OUT_do_nothing
% does not store any output, generally used during model spin-up
% S. Westermann, Jan 2021
%========================================================================


classdef OUT_do_nothing < matlab.mixin.Copyable
 

    properties
		out_index
        STRATIGRAPHY
        LATERAL
        TIMESTAMP
        MISC
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
		
	end
    
    
    methods
		
        
        function out = provide_PARA(out)         
            out.PARA.display_timestep = [];
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)

            if ~isempty(out.PARA.display_timestep)
                out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.display_timestep;
            else
                out.OUTPUT_TIME = tile.FORCING.PARA.end_time +1; %set to time that is never reached
            end
        end
        
        %---------------time integration-------------

            
        function out = store_OUT(out, tile)           
            if tile.t==out.OUTPUT_TIME 
                disp([datestr(tile.t)])
                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.display_timestep;
            end
            
        end

                %-------------param file generation-----
        function out = param_file_info(out)
            out = provide_PARA(out);

            out.PARA.STATVAR = [];
            out.PARA.options = [];
            out.PARA.class_category = 'OUT';
           
            out.PARA.default_value.display_timestep = {5};
            out.PARA.comment.display_timestep = {'timestep that model progress is displayed [days]'};

        end
        
%         function xls_out = write_excel(out)
%             XLS_OUT  Is a cell array corresponding to the class-specific content of the parameter excel file (refer to function write_controlsheet).
%             
%             xls_out = {'OUT','index',NaN,NaN;'OUT_all',1,NaN,NaN;'output_timestep',0.250000000000000,'[days]',NaN;'save_date','01.09.','provide in format dd.mm.',NaN;'save_interval',1,'[y]','if left empty, the entire output will be written out at the end';'OUT_END',NaN,NaN,NaN};
%         end
         

        
    end
end