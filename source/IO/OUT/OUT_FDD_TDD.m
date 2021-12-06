%========================================================================
% CryoGrid OUT class defining storage format of the output 
% OUT_all stores identical copies of all GROUND classses (including STATVAR, TEMP, PARA) in the
% stratigraphy for each output timestep, while lateral classes are not stored.
% The user can specify the save date and the save interval (e.g. yearly
% files), as well as the output timestep (e.g. 6 hourly). The output files
% are in Matlab (".mat") format.
% S. Westermann, T. Ingeman-Nielsen, J. Scheer, October 2020
%========================================================================


classdef OUT_FDD_TDD < matlab.mixin.Copyable
 

    properties

        TIMESTAMP
        TEMP
        PARA
        OUTPUT_TIME
        SAVE_TIME
        CONST
        STATVAR
		
	end
    
    
    methods
		
        %initialization

        
        
        function out = provide_PARA(out)         

            out.PARA.output_timestep = [];
            out.PARA.max_depth = [];
            
            out.PARA.cell_size = 0.02;
            
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)
            
			out.OUTPUT_TIME = tile.FORCING.PARA.start_time + out.PARA.output_timestep;
            out.SAVE_TIME =  tile.FORCING.PARA.end_time;

            out.STATVAR.TDD = 0;
            out.STATVAR.FDD = 0;
            out.STATVAR.time_interval = 0;
            
            out.TEMP.new_grid = [0:out.PARA.cell_size:out.PARA.max_depth]';
            out.TEMP.new_grid = (out.TEMP.new_grid(2:end,1) + out.TEMP.new_grid(1:end-1,1))./2;

        end
        
        %---------------time integration-------------
		
% 		function out = store_OUT(out, t, TOP, BOTTOM, forcing, run_number, timestep, result_path)
            
        function out = store_OUT(out, tile)           
            
             t = tile.t;
             TOP = tile.TOP; 
             BOTTOM = tile.BOTTOM;
             run_name = tile.PARA.run_name; %tile.RUN_NUMBER;
             result_path = tile.PARA.result_path;            

            
            if t>=out.OUTPUT_TIME
                disp([datestr(t)])

                T=[];
                layerThick=[];
                CURRENT = TOP.NEXT;

                while ~isequal(CURRENT, BOTTOM)
                    class_name=class(CURRENT);
                    if ~strcmp(class_name(1:4), 'SNOW')
                        T=[T; CURRENT.STATVAR.T];
                        layerThick=[layerThick; CURRENT.STATVAR.layerThick];
                    end
                    CURRENT = CURRENT.NEXT;
                end
                
                depths = [0; cumsum(layerThick)-layerThick./2];
                T_interp = interp1(depths, [T(1,1); T], out.TEMP.new_grid);
                
                out.STATVAR.TDD = out.STATVAR.TDD + T_interp .* double(T_interp>0) .* out.PARA.output_timestep;
                out.STATVAR.FDD = out.STATVAR.FDD + T_interp .* double(T_interp<0) .* out.PARA.output_timestep;
                out.STATVAR.time_interval = out.STATVAR.time_interval + out.PARA.output_timestep; 

                out.OUTPUT_TIME = out.OUTPUT_TIME + out.PARA.output_timestep;
                if t>=out.SAVE_TIME                     
                   if ~(exist([result_path run_name])==7)
                       mkdir([result_path run_name])
                   end
                   first_stable_cell = find(out.STATVAR.TDD == 0 | out.STATVAR.FDD == 0, 1);

                   out.STATVAR.TTOP = (out.STATVAR.TDD(first_stable_cell,1) + out.STATVAR.FDD(first_stable_cell,1)) ./ out.STATVAR.time_interval;
                   out.STATVAR.TTOP_depth = out.TEMP.new_grid(first_stable_cell, 1);
                   
                   save([result_path run_name '/' run_name '_OUT_FDD_TDD.mat'], 'out')
                end
            end
        end


    end
end