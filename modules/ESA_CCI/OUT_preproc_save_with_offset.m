
classdef OUT_preproc_save_with_offset < matlab.mixin.Copyable
 

    properties

        PARA
        OUTPUT_TIME
        SAVE_TIME
        START_SAVE_INTERVAL
        TEMP
        CONST
		
	end
    
    
    methods
		

        function out = provide_PARA(out)         
            out.PARA.save_date = [];
            out.PARA.save_interval = [];
            out.PARA.output_path = [];
            out.PARA.offset_in_years = [];
            
        end
        
        function out = provide_CONST(out)

        end
        
        function out = provide_STATVAR(out)

        end
		

		
		function out = finalize_init(out, tile)

 			forcing = tile.FORCING;
%             
% 			out.OUTPUT_TIME = forcing.PARA.start_time + out.PARA.output_timestep;
            if isempty(out.PARA.save_interval) || isnan(out.PARA.save_interval) 
                out.SAVE_TIME = forcing.PARA.end_time;
            else
                out.SAVE_TIME = min(forcing.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(forcing.PARA.start_time,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
            end
            out.START_SAVE_INTERVAL = forcing.PARA.start_time;
            
            
            mkdir([out.PARA.output_path tile.RUN_INFO.PARA.run_name]);
              
            out.TEMP.ERA_snowfall_downscaled =[];
            out.TEMP.ERA_melt_T_index =[];
            out.TEMP.final_av_T =[];
            out.TEMP.ERA_T_downscaled = [];
            out.TEMP.final_MODIS_weight =[];
            
        end
        
        %---------------time integration-------------
		
            
        function out = store_OUT(out, tile)          
            
            out.TEMP.ERA_snowfall_downscaled =[out.TEMP.ERA_snowfall_downscaled tile.PREPROC_CLASS.STATVAR.ERA_snowfall_downscaled];
            out.TEMP.ERA_melt_T_index =[out.TEMP.ERA_melt_T_index tile.PREPROC_CLASS.STATVAR.ERA_melt_T_index];
            out.TEMP.final_av_T =[out.TEMP.final_av_T tile.PREPROC_CLASS.STATVAR.final_av_T];
            out.TEMP.ERA_T_downscaled = [out.TEMP.ERA_T_downscaled tile.PREPROC_CLASS.STATVAR.ERA_T_downscaled];
            out.TEMP.final_MODIS_weight =[out.TEMP.final_MODIS_weight tile.PREPROC_CLASS.STATVAR.final_MODIS_weight];
            
            if tile.t  >= out.SAVE_TIME
                
                disp('saving')
                
                start_save_string = datestr(out.START_SAVE_INTERVAL, 'yyyymmdd');
                start_save_string(1:4) = num2str(str2num(start_save_string(1:4)) + out.PARA.offset_in_years); 
                end_save_string = datestr(tile.t, 'yyyymmdd');
                end_save_string(1:4) = num2str(str2num(end_save_string(1:4)) + out.PARA.offset_in_years); 
                
                variable = 'ERA_snowfall_downscaled';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' start_save_string '_' end_save_string '.nc'];
                nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_snowfall_downscaled,1) ,'y',size(out.TEMP.ERA_snowfall_downscaled,2)}, 'FillValue','disable');
                ncwrite(filename, variable, out.TEMP.ERA_snowfall_downscaled);

                variable = 'ERA_melt_T_index';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' start_save_string '_' end_save_string '.nc'];
                nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_melt_T_index,1) ,'y',size(out.TEMP.ERA_melt_T_index,2)}, 'FillValue','disable');
                ncwrite(filename, variable, out.TEMP.ERA_melt_T_index);

                variable = 'final_av_T';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' start_save_string '_' end_save_string '.nc'];
                nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.final_av_T,1) ,'y',size(out.TEMP.final_av_T,2)}, 'FillValue','disable');
                ncwrite(filename, variable, out.TEMP.final_av_T);
                
                variable = 'ERA_T_downscaled';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' start_save_string '_' end_save_string '.nc'];
                nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_T_downscaled,1) ,'y',size(out.TEMP.ERA_T_downscaled,2)}, 'FillValue','disable');
                ncwrite(filename, variable, out.TEMP.ERA_T_downscaled);
                
                variable = 'final_MODIS_weight';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' start_save_string '_' end_save_string '.nc'];
                nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.final_MODIS_weight,1) ,'y',size(out.TEMP.final_MODIS_weight,2)}, 'FillValue','disable');
                ncwrite(filename, variable, out.TEMP.final_MODIS_weight);
                
                out.TEMP.ERA_snowfall_downscaled =[];
                out.TEMP.ERA_melt_T_index =[];
                out.TEMP.final_av_T =[];
                out.TEMP.ERA_T_downscaled = [];
                out.TEMP.final_MODIS_weight =[];
                
                out.START_SAVE_INTERVAL = tile.t;
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.t,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                

            end

        end


    end
end