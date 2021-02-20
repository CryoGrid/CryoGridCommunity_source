

classdef FORCING_GROUND_ESA_CCI2 < matlab.mixin.Copyable
    
    properties
        DATA            % forcing data time series
        TEMP            % forcing data interpolated to a timestep
        PARA            % parameters
        CONST
    end
    
    
    methods
        
        
        function forcing = provide_PARA(forcing)         

            forcing.PARA.forcing_ground_folder = [];   %filename of Matlab file containing forcing data            
            forcing.PARA.start_time = []; % start time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.end_time = [];   % end time of the simulations (must be within the range of data in forcing file)
            forcing.PARA.start_overlap = [];
            forcing.PARA.end_overlap = [];
            
            forcing.PARA.preprocessed = [];

        end
        
        
        
        function forcing = provide_CONST(forcing)
            forcing.CONST.day_sec = [];
        end
        
        function forcing = provide_STATVAR(forcing)
            
        end

        
        function forcing = finalize_init(forcing, tile)
            if ~forcing.PARA.preprocessed
                
                forcing.DATA.timestamp = ncread([forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' ...
                    tile.RUN_INFO.PARA.run_name '_1.nc'], 'timestamp');
                forcing.DATA.timestamp = forcing.DATA.timestamp';
                
                variables = {'ERA_melt_bare'; 'ERA_melt_forest'; 'ERA_snowfall_downscaled'; 'ERA_T_downscaled'; 'final_av_T'; 'final_MODIS_weight'};
                
                segment_size = ceil(size(tile.RUN_INFO.STATVAR.key,1) ./ tile.RUN_INFO.PARA.num_ranks);
                start_file_number = ceil(tile.PARA.range(1) ./ segment_size);
                end_file_number = ceil(tile.PARA.range(end) ./ segment_size);
                start_pos = mod(tile.PARA.range(1), segment_size);
                start_pos = start_pos + double(start_pos==0) .* segment_size;
                end_pos = mod(tile.PARA.range(end), segment_size);
                end_pos = end_pos + double(end_pos==0) .* segment_size;
                [segment_size start_file_number end_file_number start_pos end_pos]
                
                for i=1:size(variables,1)
                    if start_file_number==end_file_number
                        file_name = [forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' num2str(start_file_number) '.nc'];
                        forcing.DATA.(variables{i,1}) = ncread(file_name, variables{i,1}, [start_pos 1], [end_pos-start_pos+1 Inf], [1 1]);
                    else
                        file_name = [forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' num2str(start_file_number) '.nc'];
                        forcing.DATA.(variables{i,1}) = ncread(file_name, variables{i,1}, [start_pos 1], [Inf Inf], [1 1]);
                        for jj=start_file_number+1:end_file_number-1
                            file_name = [forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' num2str(jj) '.nc'];
                            forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}); ncread(file_name, variables{i,1}, [1 1], [Inf Inf], [1 1])];
                        end
                        file_name = [forcing.PARA.forcing_ground_folder tile.RUN_INFO.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_' num2str(end_file_number) '.nc'];
                        forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}); ncread(file_name, variables{i,1}, [1 1], [end_pos Inf], [1 1])];
                    end
                end
                
                
                %all loaded
                %statistical downscaling
                year_list = str2num(datestr(forcing.DATA.timestamp, 'yyyy'));
                start_index = find(year_list(:,1)==forcing.PARA.start_overlap, 1);
                end_index = find(year_list(:,1)==forcing.PARA.end_overlap, 1);
                
                sum_xy = zeros(size(tile.PARA.range,1), 46); %zeros(tile.RUN_INFO.PARA.number_of_cells_per_tile,46);  %here!!! size(tile.PARA.range,1)
                sum_xx = zeros(size(tile.PARA.range,1), 46); %zeros(tile.RUN_INFO.PARA.number_of_cells_per_tile,46);
                sum_x = zeros(size(tile.PARA.range,1), 46); %zeros(tile.RUN_INFO.PARA.number_of_cells_per_tile,46);
                sum_y = zeros(size(tile.PARA.range,1), 46); %zeros(tile.RUN_INFO.PARA.number_of_cells_per_tile,46);
                n = zeros(size(tile.PARA.range,1), 46); %zeros(tile.RUN_INFO.PARA.number_of_cells_per_tile,46);
                fit_window = 5; % 5*8 = 40 days
                
                for ii=start_index:46:end_index
                    for i=1:46
                        for jj=1:fit_window  % fit routine ERA vs T_final here
                            j = mod(i+jj - 4,46)+1;
                            sum_xy(:,j) = sum_xy(:,j) + forcing.DATA.final_av_T(:,ii+i-1) .* forcing.DATA.ERA_T_downscaled(:,ii+i-1);
                            sum_xx(:,j) = sum_xx(:,j) + forcing.DATA.ERA_T_downscaled(:,ii+i-1).^2;
                            sum_x(:,j) = sum_x(:,j) + forcing.DATA.ERA_T_downscaled(:,ii+i-1);
                            sum_y(:,j) = sum_y(:,j) + forcing.DATA.final_av_T(:,ii+i-1);
                            n(:,j) =  n(:,j) + 1;
                        end
                    end
                end
                
                slope = (sum_xy - sum_x .* sum_y ./ n) ./ ( sum_xx - sum_x.^2 ./ n);
                intercept = sum_y ./ n - slope .* sum_x ./ n;
                
                forcing.DATA.slope = slope;
                forcing.DATA.intercept = intercept;
                
                for ii=1:46:start_index-46
                    for i=1:46
                        forcing.DATA.final_av_T(:,ii+i-1) = intercept(:,i) + slope(:,i) .* forcing.DATA.ERA_T_downscaled(:,ii+i-1);
                    end
                end
                
                tile.RUN_INFO.PPROVIDER.STORAGE.FORCING_GROUND_ESA_CCI.DATA = forcing.DATA;
            else
                forcing.DATA = tile.RUN_INFO.PPROVIDER.STORAGE.FORCING_GROUND_ESA_CCI.DATA;
            end
            
            
            
            forcing.PARA.start_time = datenum(forcing.PARA.start_time(1,1), forcing.PARA.start_time(2,1), forcing.PARA.start_time(3,1));
            forcing.PARA.end_time = datenum(forcing.PARA.end_time(1,1), forcing.PARA.end_time(2,1), forcing.PARA.end_time(3,1));
            
            %add last value, otherwise it chrashes when reaching the very
            %last timestamp
            forcing.DATA.ERA_melt_bare = [forcing.DATA.ERA_melt_bare forcing.DATA.ERA_melt_bare(:,end)];
            forcing.DATA.ERA_melt_forest = [forcing.DATA.ERA_melt_forest forcing.DATA.ERA_melt_forest(:,end)];            
            forcing.DATA.ERA_snowfall_downscaled = [forcing.DATA.ERA_snowfall_downscaled forcing.DATA.ERA_snowfall_downscaled(:,end)];
            forcing.DATA.final_av_T = [forcing.DATA.final_av_T forcing.DATA.final_av_T(:,end)];
            forcing.DATA.final_MODIS_weight = [forcing.DATA.final_MODIS_weight forcing.DATA.final_MODIS_weight(:,end)];
            forcing.DATA.ERA_T_downscaled = [forcing.DATA.ERA_T_downscaled forcing.DATA.ERA_T_downscaled(:,end)];
            forcing.DATA.timestamp = [forcing.DATA.timestamp forcing.PARA.end_time+1];
            
            forcing.TEMP.index = 1;
            forcing.TEMP.fraction = 0;
            while forcing.DATA.timestamp(1, forcing.TEMP.index) <= forcing.PARA.start_time
                forcing.TEMP.index = forcing.TEMP.index + 1;
            end
            forcing.TEMP.index = max(1, forcing.TEMP.index - 1);
            forcing.TEMP.number_of_substeps  =  round((forcing.DATA.timestamp(1, forcing.TEMP.index+1) - forcing.DATA.timestamp(1, forcing.TEMP.index)) .* forcing.CONST.day_sec ./ tile.timestep) ;
            forcing.TEMP.fraction = round((forcing.PARA.start_time - forcing.DATA.timestamp(1, forcing.TEMP.index)) ./ ...
                (forcing.DATA.timestamp(1, forcing.TEMP.index+1) - forcing.DATA.timestamp(1, forcing.TEMP.index)) .* forcing.TEMP.number_of_substeps);
            
            
        end
        
        function forcing = interpolate_forcing(forcing, tile)

            forcing.TEMP.surfT = (forcing.DATA.final_av_T(:,forcing.TEMP.index) + forcing.TEMP.fraction./forcing.TEMP.number_of_substeps .* ...
                (forcing.DATA.final_av_T(:, forcing.TEMP.index+1) - forcing.DATA.final_av_T(:, forcing.TEMP.index)))';
            forcing.TEMP.snowfall = (forcing.DATA.ERA_snowfall_downscaled(:,forcing.TEMP.index) + forcing.TEMP.fraction./forcing.TEMP.number_of_substeps .* ...
                (forcing.DATA.ERA_snowfall_downscaled(:, forcing.TEMP.index+1) - forcing.DATA.ERA_snowfall_downscaled(:, forcing.TEMP.index)))';
            forcing.TEMP.melt_bare = (forcing.DATA.ERA_melt_bare(:,forcing.TEMP.index) + forcing.TEMP.fraction./forcing.TEMP.number_of_substeps .* ...
                (forcing.DATA.ERA_melt_bare(:, forcing.TEMP.index+1) - forcing.DATA.ERA_melt_bare(:, forcing.TEMP.index)))';
            forcing.TEMP.melt_forest = (forcing.DATA.ERA_melt_forest(:,forcing.TEMP.index) + forcing.TEMP.fraction./forcing.TEMP.number_of_substeps .* ...
                (forcing.DATA.ERA_melt_forest(:, forcing.TEMP.index+1) - forcing.DATA.ERA_melt_forest(:, forcing.TEMP.index)))';
            
            forcing.TEMP.fraction =  forcing.TEMP.fraction + 1;
            
            if forcing.TEMP.fraction == forcing.TEMP.number_of_substeps

                forcing.TEMP.fraction = 0;
                forcing.TEMP.index = forcing.TEMP.index + 1;
                forcing.TEMP.number_of_substeps  =  round((forcing.DATA.timestamp(1, forcing.TEMP.index+1) - forcing.DATA.timestamp(1, forcing.TEMP.index)) .* forcing.CONST.day_sec ./ tile.timestep) ;
                disp(datestr(forcing.DATA.timestamp(1, forcing.TEMP.index)))
            end
            
%             posit=floor((t-forcing.DATA.timeForcing(1,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1)))+1;
%             
%             forcing.TEMP.snowfall=forcing.DATA.snowfall(posit,1)+(forcing.DATA.snowfall(posit+1,1)-forcing.DATA.snowfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.rainfall=forcing.DATA.rainfall(posit,1)+(forcing.DATA.rainfall(posit+1,1)-forcing.DATA.rainfall(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.Lin=forcing.DATA.Lin(posit,1)+(forcing.DATA.Lin(posit+1,1)-forcing.DATA.Lin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.Sin=forcing.DATA.Sin(posit,1)+(forcing.DATA.Sin(posit+1,1)-forcing.DATA.Sin(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.Tair=forcing.DATA.Tair(posit,1)+(forcing.DATA.Tair(posit+1,1)-forcing.DATA.Tair(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.wind=forcing.DATA.wind(posit,1)+(forcing.DATA.wind(posit+1,1)-forcing.DATA.wind(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.q=forcing.DATA.q(posit,1)+(forcing.DATA.q(posit+1,1)-forcing.DATA.q(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.p=forcing.DATA.p(posit,1)+(forcing.DATA.p(posit+1,1)-forcing.DATA.p(posit,1)).*(t-forcing.DATA.timeForcing(posit,1))./(forcing.DATA.timeForcing(2,1)-forcing.DATA.timeForcing(1,1));
%             
%             forcing.TEMP.rainfall = forcing.TEMP.rainfall + double(forcing.TEMP.Tair > 2) .* forcing.TEMP.snowfall;  %reassign unphysical snowfall
%             forcing.TEMP.snowfall = double(forcing.TEMP.Tair <= 2) .* forcing.TEMP.snowfall;
%             forcing.TEMP.t = t;
        end

                
    end
end