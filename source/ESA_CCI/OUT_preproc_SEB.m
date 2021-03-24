
classdef OUT_preproc_SEB < matlab.mixin.Copyable
 

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
            
            if ~(exist([out.PARA.output_path tile.RUN_INFO.PARA.run_name])==7)
                mkdir([out.PARA.output_path tile.RUN_INFO.PARA.run_name]);
            end
            %mkdir([out.PARA.output_path tile.RUN_INFO.PARA.run_name]);
              
            out.TEMP.ERA_snowfall_downscaled =[];
            out.TEMP.ERA_melt_bare =[];
            out.TEMP.ERA_melt_forest =[];
            out.TEMP.final_av_T =[];
            out.TEMP.ERA_T_downscaled = [];
            out.TEMP.final_MODIS_weight =[];
            
            out.TEMP.albedo_forest = [];
            out.TEMP.albedo_bare = [];
            
            out.TEMP.Lin = [];
            out.TEMP.Sin = [];
        end
        
        %---------------time integration-------------
		
            
        function out = store_OUT(out, tile)          
            
            out.TEMP.ERA_snowfall_downscaled =[out.TEMP.ERA_snowfall_downscaled tile.PREPROC_CLASS.STATVAR.ERA_snowfall_downscaled];
            out.TEMP.ERA_melt_bare =[out.TEMP.ERA_melt_bare tile.PREPROC_CLASS.STATVAR.ERA_melt_bare];
            out.TEMP.ERA_melt_forest =[out.TEMP.ERA_melt_forest tile.PREPROC_CLASS.STATVAR.ERA_melt_forest];
            
            out.TEMP.final_av_T =[out.TEMP.final_av_T tile.PREPROC_CLASS.STATVAR.final_av_T];
            out.TEMP.ERA_T_downscaled = [out.TEMP.ERA_T_downscaled tile.PREPROC_CLASS.STATVAR.ERA_T_downscaled];
            out.TEMP.final_MODIS_weight =[out.TEMP.final_MODIS_weight tile.PREPROC_CLASS.STATVAR.final_MODIS_weight];
            
            out.TEMP.albedo_forest = [out.TEMP.albedo_forest tile.PREPROC_CLASS.STATVAR.albedo_forest];
            out.TEMP.albedo_bare = [out.TEMP.albedo_bare tile.PREPROC_CLASS.STATVAR.albedo_bare];
            
%             out.TEMP.Lin = [out.TEMP.Lin tile.FORCING.TEMP.ERA_Sin_downscaled(:,(1:end-4))];
%             out.TEMP.Sin = [out.TEMP.Sin tile.FORCING.TEMP.ERA_Lin_downscaled(:,(1:end-4))];
            
            %if tile.t  >= out.SAVE_TIME
                
                disp('saving')

                variable = 'ERA_snowfall_downscaled';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_snowfall_downscaled,1) ,'y',size(out.TEMP.ERA_snowfall_downscaled,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.ERA_snowfall_downscaled);

                variable = 'ERA_melt_bare';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_melt_bare,1) ,'y',size(out.TEMP.ERA_melt_bare,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.ERA_melt_bare);
                
                variable = 'ERA_melt_forest';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_melt_forest,1) ,'y',size(out.TEMP.ERA_melt_forest,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.ERA_melt_forest);

                variable = 'final_av_T';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.final_av_T,1) ,'y',size(out.TEMP.final_av_T,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.final_av_T);
                
                variable = 'ERA_T_downscaled';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.ERA_T_downscaled,1) ,'y',size(out.TEMP.ERA_T_downscaled,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.ERA_T_downscaled);
                
                variable = 'final_MODIS_weight';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.final_MODIS_weight,1) ,'y',size(out.TEMP.final_MODIS_weight,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.final_MODIS_weight);
                
                variable = 'albedo_forest';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.albedo_forest,1) ,'y',size(out.TEMP.albedo_forest,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.albedo_forest);
                
                variable = 'albedo_bare';
                filename=[out.PARA.output_path tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' variable '_' datestr(out.START_SAVE_INTERVAL, 'yyyymmdd') '.nc'];
                if ~(exist(filename)==2)
                    nccreate(filename, variable, 'Dimensions', {'x',size(out.TEMP.albedo_bare,1) ,'y',size(out.TEMP.albedo_bare,2)}, 'FillValue','disable');
                end
                ncwrite(filename, variable, out.TEMP.albedo_bare);
                
%                 Lin = out.TEMP.Lin;
%                 Sin = out.TEMP.Sin;
%                 save('test2016.mat', 'Lin', 'Sin')
                
                out.TEMP.ERA_snowfall_downscaled =[];
                out.TEMP.ERA_melt_bare =[];
                out.TEMP.ERA_melt_forest =[];
                out.TEMP.final_av_T =[];
                out.TEMP.ERA_T_downscaled = [];
                out.TEMP.final_MODIS_weight =[];
                
                out.TEMP.albedo_forest = [];
                out.TEMP.albedo_bare = [];
                
                out.START_SAVE_INTERVAL = tile.t;
                out.SAVE_TIME = min(tile.FORCING.PARA.end_time,  datenum([out.PARA.save_date num2str(str2num(datestr(tile.t,'yyyy')) + out.PARA.save_interval)], 'dd.mm.yyyy'));
                
           % end

        end


    end
end