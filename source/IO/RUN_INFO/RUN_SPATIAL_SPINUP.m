%========================================================================
% CryoGrid RUN_INFO class RUN_SPATIAL_SPINUP
% RUN_INFO class for spatially distributed runs (using an appropriate 
% SPATIAL_REFERENCE class, DATA_PROVIDER classes and FORCING classes)
% which can run several TILE classes per point sequentially for model spin-up 
%
% S. westermann, Dec 2022
%========================================================================

classdef RUN_SPATIAL_SPINUP < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        SPATIAL
        TILE
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)

            run_info.PARA.number_of_cores = [];

            %run_info.PARA.number_of_cells_per_tile = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = []; %vector
            
            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
                        
        end
        
        function run_info = provide_CONST(run_info)

        end
        
        function run_info = provide_STATVAR(run_info)

        end
        
        
        function run_info = finalize_init(run_info)
            
            disp('get spatial data')
            run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
            run_info.SPATIAL.RUN_INFO = run_info;
            run_info.SPATIAL = finalize_init(run_info.SPATIAL);
        end
        
        
        function [run_info, tile] = run_model(run_info)
            
            tile = 0;
            if run_info.PARA.number_of_cores > 1
                parpool(run_info.PARA.number_of_cores)
                spmd
                    for run_number = 1:size(run_info.SPATIAL.STATVAR.key,1) %this is still wrong, does not distribute the load over workers, must be taken from ESA_CCI
                        
                        for ai=1:size(run_info.SPATIAL.ACTION,1)
                            run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, run_number); %writes the provider class
                        end
                        
                        disp(['running grid cell ' num2str(run_number)])
                        %as normal 1D run
                        for i=1:size(run_info.PARA.tile_class,1)
                            disp(['running tile number ' num2str(i)])
                            for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                                disp(['running round ' num2str(j)])
                                
                                new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                                new_tile.RUN_INFO = run_info;
                                new_tile = finalize_init(new_tile);
                                tile = new_tile;
                                run_info.TILE = tile;
                                
                                [tile.PARA.latitude tile.PARA.longitude]
                                
                                tile = run_model(tile);  %time integration
                            end
                        end
                    end
                end
                delete(gcp('nocreate'));
            else
                for run_number = 1:size(run_info.SPATIAL.STATVAR.key,1)    
                    for i=1:size(run_info.SPATIAL.ACTION,1)
                        run_info.SPATIAL.ACTION{i,1} = assign_tile_properties(run_info.SPATIAL.ACTION{i,1}, run_number); %writes the provider class
                    end
                    
                    disp(['running grid cell ' num2str(run_number)])
                    %as normal 1D run
                    for i=1:size(run_info.PARA.tile_class,1) 
                        disp(['running tile number ' num2str(i)])
                        for j=1:run_info.PARA.number_of_runs_per_tile(i,1)
                            disp(['running round ' num2str(j)])
                            
                            new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                            new_tile.RUN_INFO = run_info;
                            new_tile = finalize_init(new_tile);
                            tile = new_tile;
                            run_info.TILE = tile;
                            
                            [tile.PARA.latitude tile.PARA.longitude]
                            
                            tile = run_model(tile);  %time integration
                        end
                    end
                end
            end
            
        end
        
        
        
        %-------------param file generation-----
        function run_info = param_file_info(run_info)
            run_info = provide_PARA(run_info);

            run_info.PARA.STATVAR = [];
            run_info.PARA.class_category = 'RUN_INFO';
            run_info.PARA.default_value = [];
            run_info.PARA.comment = [];
            
            run_info.PARA.comment.number_of_cores = {'number of cores to be used for calculation'};
            run_info.PARA.default_value.number_of_cores = {2};
            
            run_info.PARA.options.tile_class.name =  'H_LIST';
            run_info.PARA.options.tile_class.entries_x = {'TILE_1D_standard' 'TILE_1D_standard'};
            
            run_info.PARA.options.tile_class_index.name =  'H_LIST'; 
            run_info.PARA.options.tile_class_index.entries_x = {1 2};
            
            run_info.PARA.options.number_of_runs_per_tile.name =  'H_LIST'; % 
            run_info.PARA.options.number_of_runs_per_tile.entries_x = {1 1};
            
            run_info.PARA.comment.projection_class = {'projection class providing providing information on the locations and additinal data for each target point'};
            
        end
        
        
            
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12));
%             
%             while ~terminate_check(run_info, 'section2_done', number_of_tiles)
%                 %for run_index = 1:number_of_tiles
%                 for run_index=[[(run_info.PARA.my_rank-1).*domains_per_worker+1:number_of_tiles] [1:(run_info.PARA.my_rank-1).*domains_per_worker]]
%                     
%                     if run_slice_yes_no(run_info, 'section2_started', 'section2_done', run_index)
%                     
%                         crap = write_check(run_info, 'section2_started', run_index);
%                         
%                         disp(['running range index ' num2str(run_index)])
%                         
%                         range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
%                         
%                         for i=1:size(run_info.PARA.tile_class,1)
%                             disp(['running tile number ' num2str(i)])
%                             for j=1:run_info.PARA.number_of_runs(i,1)
%                                 disp(['running round ' num2str(j)])
%                                 
%                                 %load the next tile from the PROVIDER
%                                 tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
%                                 tile.PARA.number_of_realizations = size(range,1);
%                                 tile.PARA.range = range;
%                                 
%                                 tile.PARA.geothermal = run_info.STATVAR.geothermal(range,1);
%                                 
%                                 
%                                 tile.RUN_INFO = run_info;
%                                 
%                                 tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                                 
%                                 tile = run_model(tile);
%                             end
%                         end
%                         crap = write_check(run_info, 'section2_done', run_index);
%                     end
%                 end
%             end
            

%         end
        
        
%         function [run_info, tile] = run_postproc(run_info)
%             
%             number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
%             %domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.num_ranks - 1e-12))
%             %domains_per_worker = max(1, ceil(number_of_tiles ./ run_info.PARA.num_ranks));
%             domains_per_worker_breaks = round([0:number_of_tiles./run_info.PARA.num_ranks:number_of_tiles]');
%             
%             for i=1:size(run_info.PARA.tile_postproc_class,1)
%                 tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_postproc_class{i,1}){run_info.PARA.tile_postproc_class_index(i,1),1});
%                 
%                 %tile.PARA.run_index = [(run_info.PARA.my_rank-1).*domains_per_worker+1:min(run_info.PARA.my_rank.*domains_per_worker, number_of_tiles)]';
%                 tile.PARA.run_index = [domains_per_worker_breaks(run_info.PARA.my_rank)+1:domains_per_worker_breaks(run_info.PARA.my_rank + 1)]';
%                 if ~isempty(tile.PARA.run_index)
%                     
%                     tile.RUN_INFO = run_info;
%                     
%                     tile = finalize_init(tile); %here, tile can still access a potentially existing tile through til.RUN_INFO.TILE
%                     
%                     tile = run_model(tile);
%                 end
%             end
%        
%             %terminate parallel environment 
%             if run_info.PARA.parallelized == 1
%                 NMPI_Finalize(); % End the MPI communication
%             end
%         end
 
        
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
        
%         function run_info = customize(run_info)
% 
%             
%         end
%         
%         
%          function run_info = parallelize_preproc(run_info)
%             number_of_years = [];
%             load_index = [1; 1; 1.8];
%             number_of_cores = run_info.PARA.num_ranks;
%             
%             my_core = run_info.PARA.my_rank;
%             
%             number_of_slices = 46;
%             start_year= [];
% 
%             run_info.PARA.tile_active = [0;0;0];
%             
%             number_of_slices = 46;
%             start_year= [];
%             
%             for i=1:size(run_info.PARA.tile_preproc_class_index,1)
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
%                 forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
%                 number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
%                 start_year= [start_year; forcing.PARA.start_time(1)];
%             end
%             
%             load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
%             number_of_slices_per_tile = number_of_years .* number_of_slices;
%             load_per_tile = number_of_years .* number_of_slices .* load_index;
%             load_per_tile_acc = cumsum(load_per_tile);
%             
% 
%             load_start = load_per_core .* (my_core-1)+1;
%             load_end = load_per_core .* my_core;
%             
%             tile_num = 1;
%             slice_count = 0;
%             while tile_num<=3
%                 if load_start > load_per_tile_acc(tile_num)
%                     tile_num=tile_num+1;
%                 else
%                     break
%                 end
%             end
%             tile_num_start = tile_num;
%             
%             tile_num = 1;
%             slice_count = 0;
%             while tile_num<=3
%                 if load_end > load_per_tile_acc(tile_num)
%                     tile_num=tile_num+1;
%                 else
%                     break
%                 end
%             end
%             tile_num_end = tile_num;
%             run_info.PARA.tile_active(tile_num_start:tile_num_end) = run_info.PARA.tile_active(tile_num_start:tile_num_end)+1;
%             
%             load_per_tile_acc = [0; load_per_tile_acc];
%             slice_count_start = round((load_start - load_per_tile_acc(tile_num_start)) ./ load_index(tile_num_start));
%             if slice_count_start > 1
%                 %new_start_time = datenum(start_year(tile_num_start,1) + floor(slice_count_start./46),1,1) + (mod(slice_count_start, 46)-1).*8;  
%                 new_start_time = datenum(start_year(tile_num_start,1) + floor((slice_count_start-1)./46),1,1) + (mod(slice_count_start-1, 46)).*8;
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_start,1}){run_info.PARA.tile_preproc_class_index(tile_num_start,1),1};
%                 [year,month,day,~,~,~] = datevec(new_start_time);
%                 run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day];
%                 disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
%             end
%             slice_count_end =  round((load_end - load_per_tile_acc(tile_num_end)) ./ load_index(tile_num_end))+2; %add the 2 for overlap
%             if slice_count_end < number_of_slices_per_tile(tile_num_end)
%                 %new_end_time = datenum(start_year(tile_num_end,1) + floor(slice_count_end./46),1,1) + (mod(slice_count_end, 46)-1).*8;
%                 new_end_time = datenum(start_year(tile_num_end,1) + floor((slice_count_end-1)./46),1,1) + (mod(slice_count_end-1, 46)).*8;
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num_end,1}){run_info.PARA.tile_preproc_class_index(tile_num_end,1),1};
%                 [year,month,day,~,~,~] = datevec(new_end_time);
%                 run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day];
%                 disp([num2str(run_info.PARA.my_rank) ' ' num2str(year) ' ' num2str(month) ' ' num2str(day)])
%             end
% 
%          end
%         
%         
%         
%         function run_info = parallelize_preproc_old(run_info)
%             number_of_years = [];
%             load_index = [1; 1; 1.8];
%             number_of_cores = run_info.PARA.num_ranks;
%             
%             my_core = run_info.PARA.my_rank;
%             
%             number_of_slices = 46;
%             start_year= [];
%             
% 
%             run_info.PARA.tile_active = [0;0;0];
%             
%             number_of_slices = 46;
%             start_year= [];
%             
%             for i=1:size(run_info.PARA.tile_preproc_class_index,1)
%                 tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{i,1}){run_info.PARA.tile_preproc_class_index(i,1),1};
%                 forcing = run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1};
%                 number_of_years = [number_of_years; forcing.PARA.end_time(1)-forcing.PARA.start_time(1)+1];
%                 start_year= [start_year; forcing.PARA.start_time(1)];
%             end
%             
%             load_per_core = sum(number_of_years .* number_of_slices .* load_index) ./ number_of_cores;
%             number_of_slices_per_tile = number_of_years .* number_of_slices;
%             load_per_tile = number_of_years .* number_of_slices .* load_index;
%             load_per_tile_acc = cumsum(load_per_tile);
%             
%             result =[];
%             tile_number = 1;
%             core = 1;
%             slice_count = 0;
%             end_last_slice = 0;
%             load_count = 0;
%             
%             
%             while core <= number_of_cores %tile_number<=3
%                 load_count = load_count + load_per_core;
%                 if load_count < load_per_tile_acc(tile_number)
%                     slice_count = slice_count + load_per_core./load_index(tile_number);
%                     result=[result; [core tile_number round(slice_count)]];
%                 else
%                     while tile_number<3
%                         remaining_load = load_count - load_per_tile_acc(tile_number);
%                         if  remaining_load >0
%                             result=[result; [core tile_number number_of_slices_per_tile(tile_number)]];
%                         else
%                             break
%                         end
%                         tile_number = tile_number+1;
%                         slice_count =  remaining_load./load_index(tile_number);
%                         if round(slice_count) <= number_of_slices_per_tile(tile_number)
%                             result=[result; [core tile_number round(slice_count)]];
%                         end
%                     end
%                 end
%                 
%                 core = core +1;
%             end
%             
%             result2=[];
%             ti=0;
%             for i=1:size(result,1)
%                 
%                 if ti ~= result(i,2)
%                     result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%                         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%                         datenum(start_year(result(i,2)),1,1)]];
%                     ti = result(i,2);
%                 else
%                     result2 = [result2; [start_year(result(i,2)) + floor(result(i,3)./46) mod(result(i,3),46).*8 + 1 ...
%                         datenum(start_year(result(i,2)) + floor(result(i,3)./46),1,1) + mod(result(i,3),46).*8-1  ...
%                         datenum(start_year(result(i-1,2)) + floor(result(i-1,3)./46),1,1) + mod(result(i-1,3),46).*8]];
%                 end
%             end
%             % result2 = [result2 [datenum(start_year(1,1), 1, 1); result2(1:end-1,3)]];
%             % result2(:,3) = result2(:,3)-1;
%             
%             for i=1:size(result,1)
%                 if my_core == result(i,1)
%                     tile_num =result(i,2);
%                     run_info.PARA.tile_active(tile_num,1) = 1;
% 
%                     tile = run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_preproc_class{tile_num,1}){run_info.PARA.tile_preproc_class_index(tile_num,1),1};
%                     [year,month,day,~,~,~] = datevec(result2(i,4));
%                     run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.start_time = [year; month; day]; 
%                     disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
%                     [year,month,day,~,~,~] = datevec(result2(i,3));
%                     run_info.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1}.PARA.end_time = [year; month; day]; 
%                     disp([num2str(run_info.PARA.my_rank) num2str(year) ' ' num2str(month) ' ' num2str(day)])
% %                     new_start_time = datestr(result2(i,4))
% %                     new_end_time = datestr(result2(i,3))
%                 end
%                 
%             end
%             
%         end
        
        
    end
end



