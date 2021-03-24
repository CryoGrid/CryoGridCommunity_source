
classdef RUN_ESA_CCI < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.run_name = [];
            run_info.PARA.parallelized = [];
            
            run_info.PARA.checkfolder = [];
            
            run_info.PARA.tile_preproc_class = [];
            run_info.PARA.tile_preproc_class_index = [];
            
            run_info.PARA.tile_rearrange_files_class = [];
            run_info.PARA.tile_rearrange_files_class_index = [];

            run_info.PARA.number_of_cells_per_tile = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs = [];
            
            run_info.PARA.tile_postproc_class = [];
            run_info.PARA.tile_postproc_class_index = [];
            
            run_info.PARA.base_projection_class = [];
            run_info.PARA.base_projection_class_index = [];
            
            run_info.PARA.mask_class = []; %generally several classes - array
            run_info.PARA.mask_class_index = [];
            
            run_info.PARA.DEM_class = [];
            run_info.PARA.DEM_class_index = [];
            
            run_info.PARA.landcover_class = [];
            run_info.PARA.landcover_class_index = [];
            
            run_info.PARA.geothermal_class = [];
            run_info.PARA.geothermal_class_index = [];
            
            
            

            
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)


            if run_info.PARA.parallelized == 1
                addpath('/cluster/software/MATLAB/2019a/NMPI/version13_intel2019a/');
                NMPI_Init(); % Init the MPI communication between MPI processes
                run_info.PARA.num_ranks =  NMPI_Comm_size(); % Get number of MPI processes.%
                run_info.PARA.my_rank   =  NMPI_Comm_rank()+1; % Get my MPI process ID (from 0 to num_ranks-1)
                if ~(exist([run_info.PARA.checkfolder run_info.PARA.run_name])==7)
                    mkdir([run_info.PARA.checkfolder run_info.PARA.run_name])
                end
                run_info.PARA.checkfolder = [run_info.PARA.checkfolder run_info.PARA.run_name '/'];
            else
                run_info.PARA.num_ranks =  1;
                run_info.PARA.my_rank   =  1;
                if ~(exist([run_info.PARA.checkfolder run_info.PARA.run_name])==7)
                    mkdir([run_info.PARA.checkfolder run_info.PARA.run_name])
                end
                run_info.PARA.checkfolder = [run_info.PARA.checkfolder run_info.PARA.run_name '/'];
            end

            disp('get spatial reference')
            %provided by coordinate system MODIS LST classes
            run_info.STATVAR.spatial = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.base_projection_class){run_info.PARA.base_projection_class_index,1});
            run_info.STATVAR.spatial.RUN_INFO = run_info;
            run_info.STATVAR.spatial = finalize_init(run_info.STATVAR.spatial);
            
            run_info.STATVAR.latitude = run_info.STATVAR.spatial.STATVAR.latitude;
            run_info.STATVAR.longitude = run_info.STATVAR.spatial.STATVAR.longitude;
            run_info.STATVAR.key = run_info.STATVAR.spatial.STATVAR.key;
            run_info.STATVAR.list_of_MODIS_tiles = run_info.STATVAR.spatial.STATVAR.list_of_MODIS_tiles; %use to loop over several tiles
            run_info.PARA.total_number_of_cells = size(run_info.STATVAR.key,1);
            run_info.STATVAR.properties = run_info.STATVAR.spatial.STATVAR.properties;
            
            disp('get altitudes')
            run_info.PARA.dem = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.DEM_class){run_info.PARA.DEM_class_index,1});
            run_info.PARA.dem = finalize_init(run_info.PARA.dem);
            run_info.STATVAR.altitude = get_altitudes(run_info.PARA.dem, run_info);

            disp('get geothermal')
            run_info.PARA.geothermal = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.geothermal_class){run_info.PARA.geothermal_class_index,1});
            run_info.PARA.geothermal = finalize_init(run_info.PARA.geothermal);
            run_info.STATVAR.geothermal = get_geothermal_heatflux(run_info.PARA.geothermal, run_info);
            
            disp('get landcover')
            run_info.PARA.landcover = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.landcover_class){run_info.PARA.landcover_class_index,1});
            run_info.PARA.landcover = finalize_init(run_info.PARA.landcover);
            run_info.STATVAR.landcover = get_landcover(run_info.PARA.landcover, run_info);
            
            data_package = run_info.STATVAR;
            data_package.spatial = [];
            save([run_info.PPROVIDER.PARA.result_path run_info.PPROVIDER.PARA.run_name '/' run_info.PARA.run_name '_metadata.mat'], 'data_package');


        end
        
        function [run_info, tile] = run_preproc(run_info)
            
            run_info = parallelize_preproc(run_info);
           

            for i=1:size(run_info.PARA.tile_preproc_class,1)
                if run_info.PARA.tile_active(i,1) ==1
                    tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1});
                    tile.RUN_INFO = run_info;
                    tile = finalize_init(tile);
                    tile = run_model(tile);
                end
            end
            success = write_check(run_info, 'section1_done', run_info.PARA.my_rank);
            
            while ~terminate_check(run_info, 'section1_done', run_info.PARA.num_ranks)
            
            end

                %replaces %NMPI_Barrier();
                
            %rearrange forcing files so that entire time series is
            %comprised in a single nc-file
            for i=1:size(run_info.PARA.tile_rearrange_files_class, 1)
                disp(i)
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_rearrange_files_class{i,1}){run_info.PARA.tile_rearrange_files_class_index(i,1),1});
                tile.RUN_INFO = run_info;
                tile = finalize_init(tile);
                tile = run_model(tile);
                disp('done')
                success = write_check(run_info, ['rearrange_files' num2str(i) '_done'], run_info.PARA.my_rank);
                while ~terminate_check(run_info, ['rearrange_files' num2str(i) '_done'], run_info.PARA.num_ranks)
                end
            end
                
        end
        
        
        function [run_info, tile] = run_model(run_info)
            

            
            number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
            domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12));
            
            while ~terminate_check(run_info, 'section2_done', number_of_tiles)
                %for run_index = 1:number_of_tiles
                for run_index=[[(run_info.PARA.my_rank-1).*domains_per_worker+1:number_of_tiles] [1:(run_info.PARA.my_rank-1).*domains_per_worker]]
                    
                    if run_slice_yes_no(run_info, 'section2_started', 'section2_done', run_index)
                    
                        crap = write_check(run_info, 'section2_started', run_index);
                        
                        disp(['running range index ' num2str(run_index)])
                        
                        range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
                        
                        for i=1:size(run_info.PARA.tile_class,1)
                            disp(['running tile number ' num2str(i)])
                            for j=1:run_info.PARA.number_of_runs(i,1)
                                disp(['running round ' num2str(j)])
                                
                                %load the next tile from the PROVIDER
                                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                                tile.PARA.number_of_realizations = size(range,1);
                                tile.PARA.range = range;
                                
                                tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
                                
                                
                                tile.RUN_INFO = run_info;
                                
                                tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
                                
                                tile = run_model(tile);
                            end
                        end
                        crap = write_check(run_info, 'section2_done', run_index);
                    end
                end
            end
            

        end
        
        
        function [run_info, tile] = run_postproc(run_info)
            
            number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
            %domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12))
            %domains_per_worker = max(1, ceil(number_of_tiles ./ run_info.PARA.num_ranks));
            domains_per_worker_breaks = round([0:number_of_tiles./run_info.PARA.num_ranks:number_of_tiles]');
            
            for i=1:size(run_info.PARA.tile_postproc_class,1)
                tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_postproc_class{i,1}){run_info.PARA.tile_postproc_class_index(i,1),1});
                
                %tile.PARA.run_index = [(run_info.PARA.my_rank-1).*domains_per_worker+1:min(run_info.PARA.my_rank.*domains_per_worker, number_of_tiles)]';
                tile.PARA.run_index = [domains_per_worker_breaks(run_info.PARA.my_rank)+1:domains_per_worker_breaks(run_info.PARA.my_rank + 1)]';

                tile.RUN_INFO = run_info;
                
                tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
                
                tile = run_model(tile);
            end
       
            %terminate parallel environment 
            if run_info.PARA.parallelized == 1
                NMPI_Finalize(); % End the MPI communication
            end
        end
 
        
%         
%         function [run_info, tile] = run_model(run_info)
%             
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12));
%             
%             %parallelize this
%             for run_index = 1:number_of_tiles
%                 
%                 disp(['running range index ' num2str(run_index)])
%                 
%                 range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
%                 
%                 for i=1:size(run_info.PARA.tile_class,1)
%                     disp(['running tile number ' num2str(i)])
%                     for j=1:run_info.PARA.number_of_runs(i,1)
%                         disp(['running round ' num2str(j)])
%                         
%                         %load the next tile from the PROVIDER
%                         tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
%                         tile.PARA.number_of_realizations = size(range,1);
%                         tile.PARA.range = range;
%                         
%                         tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
%                         
%                         
%                         %REMOVE
% %                         [~, pos] = max(run_info.STATVAR.landcover(range,1),[], 2);
% %                         tile.PARA.stratigraphy = pos(1,1) .*0 +1;
%                         %REMOVE
%                         
%                         %tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class){run_info.PARA.tile_class_index,1});
%                         tile.RUN_INFO = run_info;
%                         
%                         tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                        
%                         tile = run_model(tile);
%                     end
%                 end
%             end
%             
%             if run_info.PARA.parallelized == 1
%                 NMPI_Finalize(); % End the MPI communication
%             end
%         end
 
        
        
        %non-mandatory
        
        function run_info = customize(run_info)

            
        end
        
        
         function run_info = parallelize_preproc(run_info)
            number_of_years = [];
            load_index = [1; 1; 1.8];
            number_of_cores = run_info.PARA.num_ranks;
            
            my_core = run_info.PARA.my_rank;
            
            number_of_slices = 46;
            start_year= [];

            run_info.PARA.tile_active = [0;0;0];
            
            number_of_slices = 46;
            start_year= [];
            
            for i=1:size(run_info.PARA.tile_preproc_class_index,1)
                tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
                forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
                number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
                start_year= [start_year; forcing.PARA.start_time(1)];
            end
            
            load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
            number_of_slices_per_tile = number_of_years .* number_of_slices;
            load_per_tile = number_of_years .* number_of_slices .* load_index;
            load_per_tile_acc = cumsum(load_per_tile);
            

            load_start = load_per_core .* (my_core-1)+1;
            load_end = load_per_core .* my_core;
            
            tile_num = 1;
            slice_count = 0;
            while tile_num<=3
                if load_start > load_per_tile_acc(tile_num)
                    tile_num=tile_num+1;
                else
                    break
                end
            end
            tile_num_start = tile_num;
            
            tile_num = 1;
            slice_count = 0;
            while tile_num<=3
                if load_end > load_per_tile_acc(tile_num)
                    tile_num=tile_num+1;
                else
                    break
                end
            end
            tile_num_end = tile_num;
            run_info.PARA.tile_active(tile_num_start:tile_num_end) = run_info.PARA.tile_active(tile_num_start:tile_num_end)+1;
            
            load_per_tile_acc = [0; load_per_tile_acc];
            slice_count_start = round((load_start - load_per_tile_acc(tile_num_start)) ./ load_index(tile_num_start));
            if slice_count_start > 1
                %new_start_time = datenum(start_year(tile_num_start,1) + floor(slice_count_start./46),1,1) + (mod(slice_count_start, 46)-1).*8;  
                new_start_time = datenum(start_year(tile_num_start,1) + floor((slice_count_start-1)./46),1,1) + (mod(slice_count_start-1, 46)).*8;
                tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_start,1}){run_info.PARA.tile_preproc_class_index(tile_num_start,1),1};
                [year,month,day,~,~,~] = datevec(new_start_time);
                run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day];
                disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
            end
            slice_count_end =  round((load_end - load_per_tile_acc(tile_num_end)) ./ load_index(tile_num_end))+2; %add the 2 for overlap
            if slice_count_end < number_of_slices_per_tile(tile_num_end)
                %new_end_time = datenum(start_year(tile_num_end,1) + floor(slice_count_end./46),1,1) + (mod(slice_count_end, 46)-1).*8;
                new_end_time = datenum(start_year(tile_num_end,1) + floor((slice_count_end-1)./46),1,1) + (mod(slice_count_end-1, 46)).*8;
                tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_end,1}){run_info.PARA.tile_preproc_class_index(tile_num_end,1),1};
                [year,month,day,~,~,~] = datevec(new_end_time);
                run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day];
                disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
            end

         end
        
        
        
        function run_info = parallelize_preproc_old(run_info)
            number_of_years = [];
            load_index = [1; 1; 1.8];
            number_of_cores = run_info.PARA.num_ranks;
            
            my_core = run_info.PARA.my_rank;
            
            number_of_slices = 46;
            start_year= [];
            

            run_info.PARA.tile_active = [0;0;0];
            
            number_of_slices = 46;
            start_year= [];
            
            for i=1:size(run_info.PARA.tile_preproc_class_index,1)
                tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
                forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
                number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
                start_year= [start_year; forcing.PARA.start_time(1)];
            end
            
            load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
            number_of_slices_per_tile = number_of_years .* number_of_slices;
            load_per_tile = number_of_years .* number_of_slices .* load_index;
            load_per_tile_acc = cumsum(load_per_tile);
            
            result =[];
            tile_number = 1;
            core = 1;
            slice_count = 0;
            end_last_slice = 0;
            load_count = 0;
            
            
            while core <= number_of_cores %tile_number<=3
                load_count = load_count + load_per_core;
                if load_count < load_per_tile_acc(tile_number)
                    slice_count = slice_count + load_per_core./load_index(tile_number);
                    result=[result; [core tile_number round(slice_count)]];
                else
                    while tile_number<3
                        remaining_load = load_count - load_per_tile_acc(tile_number);
                        if  remaining_load >0
                            result=[result; [core tile_number number_of_slices_per_tile(tile_number)]];
                        else
                            break
                        end
                        tile_number = tile_number+1;
                        slice_count =  remaining_load./load_index(tile_number);
                        if round(slice_count) <= number_of_slices_per_tile(tile_number)
                            result=[result; [core tile_number round(slice_count)]];
                        end
                    end
                end
                
                core = core +1;
            end
            
            result2=[];
            ti=0;
            for i=1:size(result,1)
                
                if ti ~= result(i,2)
                    result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
                        datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
                        datenum(start_year(result(i,2)),1,1)]];
                    ti = result(i,2);
                else
                    result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
                        datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
                        datenum(start_year(result(i-1,2)) + floor(result(i-1,3)./46),1,1) + mod(result(i-1,3),46).*8]];
                end
            end
            % result2 = [result2 [datenum(start_year(1,1), 1, 1); result2(1:end-1,3)]];
            % result2(:,3) = result2(:,3)-1;
            
            for i=1:size(result,1)
                if my_core == result(i,1)
                    tile_num =result(i,2);
                    run_info.PARA.tile_active(tile_num,1) = 1;

                    tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num,1}){run_info.PARA.tile_preproc_class_index(tile_num,1),1};
                    [year,month,day,~,~,~] = datevec(result2(i,4));
                    run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day]; 
                    disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
                    [year,month,day,~,~,~] = datevec(result2(i,3));
                    run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day]; 
                    disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
%                     new_start_time = datestr(result2(i,4))
%                     new_end_time = datestr(result2(i,3))
                end
                
            end
            
        end
        
        function success = write_check(run_info, tag, n)
            dlmwrite([run_info.PARA.checkfolder tag num2str(n) '.check'], []);
            success = 1;
        end
        
        
        function success = run_slice_yes_no(run_info, start_tag, end_tag, n)
            if exist([run_info.PARA.checkfolder end_tag num2str(n) '.check'])==2
                success = 0;
            elseif ~(exist([run_info.PARA.checkfolder start_tag num2str(n) '.check'])==2)
                success = 1;
            else
                finfo = dir([run_info.PARA.checkfolder start_tag num2str(n) '.check']);
                if now - finfo.datenum < 6/24  %file could still be written by other worker
                    success = 0;
                else
                    success=1;
                end
            end
        end

        
        function success = terminate_check(run_info, end_tag, n_max)
            finfo = dir([run_info.PARA.checkfolder end_tag '*.check']);
            if size(finfo,1) >= n_max %number of check files is equal n_max
                success=1;
            else
                success = 0;
                disp(['core ' num2str(run_info.PARA.my_rank) ' pausing'])
                pause(5)
            end
        end
        
        
        
        
    end
end



