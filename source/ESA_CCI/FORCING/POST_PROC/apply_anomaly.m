%========================================================================
% CryoGrid FORCING post-processing 
%
%
% Authors:
% S. Westermann, January 2023
%
%========================================================================

classdef apply_anomaly < FORCING_base 
    
    properties
        
    end
    
    methods
        function post_proc = provide_PARA(post_proc)
            
            post_proc.PARA.anomaly_file_path = [];
            post_proc.PARA.anomaly_file = [];

            post_proc.PARA.post_proc_class = [];  
            post_proc.PARA.post_proc_class_index = [];
            post_proc.PARA.year_range = [];
            post_proc.PARA.annual = [];
            
        end
        
        
        function post_proc = provide_CONST(post_proc)
%             post_proc.CONST.L_f = []; 
%             post_proc.CONST.sigma = [];
%             post_proc.CONST.day_sec = [];
%             post_proc.CONST.Tmfw = [];
        end
        
        
        function post_proc = provide_STATVAR(post_proc)
            
        end
        
        
        function post_proc = finalize_init(post_proc, tile)

            for i=1:size(post_proc.PARA.post_proc_class,1)
                post_proc.TEMP.post_proc{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(post_proc.PARA.post_proc_class{i,1}){post_proc.PARA.post_proc_class_index(i,1),1});
                post_proc.TEMP.post_proc{i,1} = finalize_init(post_proc.TEMP.post_proc{i,1}, tile);
            end
            
        end
        
        
        function forcing = post_process(post_proc, forcing, tile)
            
            if forcing.TEMP.current_year >= post_proc.PARA.year_range(1) && forcing.TEMP.current_year <= post_proc.PARA.year_range(end)

                variables = fieldnames(forcing.TEMP); %store original data
                for i=1:size(variables, 1)
                    TEMP.(variables{i,1}) = forcing.TEMP.(variables{i,1});
                end
                
                %apply anomaly
                forcing = apply_anomaly_from_file(post_proc, forcing, tile);
                
                %run the actual processing
                for i=1:size(post_proc.TEMP.post_proc,1)
                    post_proc.TEMP.post_proc{i,1}.TEMP.offset_years = forcing.PARA.number_of_spin_up_years;
                    forcing = post_process(post_proc.TEMP.post_proc{i,1}, forcing, tile);
                end
                
                %reassign original values
                forcing.TEMP = TEMP;
                
            end
        end
        
        
        function forcing = apply_anomaly_from_file(post_proc, forcing, tile)
            
            target_lat = tile.PARA.latitude;
            target_lon = tile.PARA.longitude;
            
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'ERA_lat')
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'ERA_lon')
            
            max_target_lat = min(round(max(target_lat) /2) *2 + 2, 90) ;  %ERA5 in 0.25 degree resolution
            min_target_lat= max(round(min(target_lat)/2) *2 - 2, -90);
            
            target_lon(target_lon < 0) = target_lon(target_lon < 0) + 360; %change from -180 to +180 to 0 to 360
            max_target_lon = min(round(max(target_lon)/2) *2 + 2, 360);
            min_target_lon = max(round(min(target_lon)/2) *2 - 2, 0);
            
            lat_index_start = find(ERA_lat(:,1) == max_target_lat);
            lat_index_end = find(ERA_lat(:,1) == min_target_lat);
            lon_index_start = find(ERA_lon(:,1) == min_target_lon);
            lon_index_end = find(ERA_lon(:,1) == max_target_lon);
            
            
            [ERA_lon_mesh, ERA_lat_mesh] = meshgrid(ERA_lon(lon_index_start:lon_index_end), ERA_lat(lat_index_start:lat_index_end));
            
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'anomaly_airT_interp')
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'anomaly_precip_interp')
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'anomaly_LW_interp')
            load([post_proc.PARA.anomaly_file_path post_proc.PARA.anomaly_file], 'anomaly_SW_interp')
            
            for jj=1:size(anomaly_airT_interp,3)
                anomaly_airT_interp2=squeeze(double(anomaly_airT_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end,jj)));
                anomaly_precip_interp2=squeeze(double(anomaly_precip_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end,jj)));
                anomaly_LW_interp2=squeeze(double(anomaly_LW_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end,jj)));
                anomaly_SW_interp2=squeeze(double(anomaly_SW_interp(lon_index_start:lon_index_end, lat_index_start:lat_index_end,jj)));
                
                anomaly_temperature = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_airT_interp2', target_lon, target_lat, 'cubic');
                anomaly_precip = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_precip_interp2', target_lon, target_lat, 'cubic');
                anomaly_LW = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_LW_interp2', target_lon, target_lat, 'cubic');
                anomaly_SW = interp2(ERA_lon_mesh ,ERA_lat_mesh, anomaly_SW_interp2', target_lon, target_lat, 'cubic');
                
                averaging_period = ceil(365./size(anomaly_airT_interp,3));
                range_doy = (jj-1).*averaging_period.*4+1:min(jj*averaging_period*4, size(forcing.TEMP.ERA_T_downscaled,2));
                forcing.TEMP.ERA_T_downscaled(:, range_doy) = forcing.TEMP.ERA_T_downscaled(:, range_doy) + anomaly_temperature;
                forcing.TEMP.ERA_precip_downcaled(:, range_doy) = forcing.TEMP.ERA_precip_downcaled(:, range_doy) .* anomaly_precip;
                
                forcing.TEMP.ERA_Lin_downscaled(:, range_doy) = forcing.TEMP.ERA_Lin_downscaled(:, range_doy) + anomaly_LW;
                forcing.TEMP.ERA_Sin_downscaled(:, range_doy) = forcing.TEMP.ERA_Sin_downscaled(:, range_doy) .* anomaly_SW;
            end
        end
        
        
        
%                 %-------------param file generation-----
%         function post_proc = param_file_info(post_proc)
%             post_proc = provide_PARA(post_proc);
% 
%             post_proc.PARA.STATVAR = [];
%             post_proc.PARA.class_category = 'FORCING POST_PROCESSING';
%             post_proc.PARA.options = [];
%             
%             post_proc.PARA.eliminate_fraction = [];
%             post_proc.PARA.survive_fraction = [];
%                         
%             post_proc.PARA.default_value.window_size = {7};
%             post_proc.PARA.comment.window_size = {'window size in days within which precipitation is reallocated'};
%             
%             post_proc.PARA.default_value.eliminate_fraction = {0.5};
%             post_proc.PARA.comment.eliminate_fraction = {'fraction of smallest precipitation events (= timestamps with precipitation) that is reallocated to larger events'};
%             
%             post_proc.PARA.default_value.survive_fraction = {0.5};  
%             post_proc.PARA.comment.survive_fraction = {'fraction of largest precipitation events (= timestamps with precipitation) that the small events are reallocated to'};
%             
%         end
        
    end
    
end