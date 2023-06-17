
classdef RUN_ESA_CCI_slurmArray < matlab.mixin.Copyable
    
    properties
        PPROVIDER
        PARA
        CONST
        STATVAR
        TILE
        SPATIAL
    end
    
    
    methods
        
        
        function run_info = provide_PARA(run_info)
            
            run_info.PARA.run_name = [];
            run_info.PARA.preprocess_only = 0;
            
            run_info.PARA.checkfolder = [];
            
            run_info.PARA.number_of_cells_per_tile = [];
            
            run_info.PARA.tile_class = [];
            run_info.PARA.tile_class_index = [];
            run_info.PARA.number_of_runs_per_tile = [];
            
            run_info.PARA.projection_class = [];
            run_info.PARA.projection_class_index = [];
             

        end
        
        function run_info = provide_CONST(run_info)
            
        end
        
        function run_info = provide_STATVAR(run_info)
            
        end
        
        
        function run_info = finalize_init(run_info)
            
            run_info.PARA.parallelized = 0;
            
            if ~(exist([run_info.PARA.checkfolder run_info.PARA.run_name])==7)
                mkdir([run_info.PARA.checkfolder run_info.PARA.run_name])
            end
            
            run_info.PARA.checkfolder = [run_info.PARA.checkfolder run_info.PARA.run_name '/'];
            
            disp('get spatial data')
            run_info.SPATIAL = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.projection_class){run_info.PARA.projection_class_index,1});
            run_info.SPATIAL.RUN_INFO = run_info;
            run_info.SPATIAL = finalize_init(run_info.SPATIAL);
            if run_info.PARA.worker_number == 1
                if ~(exist([run_info.PPROVIDER.PARA.result_path run_info.PARA.run_name '/'] )==7)
                    mkdir([run_info.PPROVIDER.PARA.result_path run_info.PARA.run_name '/'])
                end
                spatial_reference = copy(run_info.SPATIAL);
                spatial_reference.RUN_INFO = [];
                spatial_reference.ACTION = [];
                spatial_reference.PARA.result_path = run_info.PPROVIDER.PARA.result_path;
                spatial_reference.PARA.run_name = run_info.PARA.run_name;
                save([run_info.PPROVIDER.PARA.result_path run_info.PARA.run_name '/' run_info.PARA.run_name '.mat'], 'spatial_reference')
                spatial_reference = [];
                
            end
            
            run_info.PARA.total_number_of_cells = size(run_info.SPATIAL.STATVAR.key,1);
            
        end
        
        
        function [run_info, tile] = run_model(run_info)
            
            tile = 0;
            
            if ~run_info.PARA.preprocess_only %do nothing if this variable is set, can be used to enable preprosessing  with batch script and "wait" command
                
                number_of_tiles = ceil(run_info.PARA.total_number_of_cells ./ run_info.PARA.number_of_cells_per_tile);
                domains_per_worker = max(1, floor(number_of_tiles ./ run_info.PARA.number_of_cores - 1e-12));
                %spmd must be placed here
                
                
                for run_index=[[(run_info.PARA.worker_number-1).*domains_per_worker+1:number_of_tiles] [1:(run_info.PARA.worker_number-1).*domains_per_worker]]
                    
                    if run_slice_yes_no(run_info, 'section_started', 'section_done', run_index)
                        
                        crap = write_check(run_info, 'section_started', run_index);
                        
                        disp(['running range index ' num2str(run_index)])
                        
                        range = [(run_index-1).*run_info.PARA.number_of_cells_per_tile+1:min(run_index.*run_info.PARA.number_of_cells_per_tile, run_info.PARA.total_number_of_cells)]';
                        
                        for ai=1:size(run_info.SPATIAL.ACTION,1)
                            run_info.SPATIAL.ACTION{ai,1} = assign_tile_properties(run_info.SPATIAL.ACTION{ai,1}, range); %writes the provider class
                        end
                        
                        for i=1:size(run_info.PARA.tile_class,1)
                            disp(['running tile number ' num2str(i)])
                            for j=1:run_info.PARA.number_of_runs_per_tile (i,1)
                                disp(['running round ' num2str(j)])
                                
                                new_tile = copy(run_info.PPROVIDER.CLASSES.(run_info.PARA.tile_class{i,1}){run_info.PARA.tile_class_index(i,1),1});
                                new_tile.RUN_INFO = run_info;
                                new_tile = finalize_init(new_tile);
                                tile = new_tile;
                                run_info.TILE = tile;
                                
                                tile = run_model(tile);
                            end
                        end
                        crap = write_check(run_info, 'section_done', run_index);
                    end
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
        
        
    end
end



