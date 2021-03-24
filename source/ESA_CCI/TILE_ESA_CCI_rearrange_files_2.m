
classdef TILE_ESA_CCI_rearrange_files_2 < matlab.mixin.Copyable
    
    properties
        PARA
        RUN_INFO
        BUILDER
        FORCING
        CONST
        OUT        
        PREPROC_CLASS
       
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)
                        
            tile.PARA.file_path = [];
            tile.PARA.in_name = [];
            %tile.PARA.out_name = [];

        end
        
        function tile = provide_CONST(tile)

        end
        
        function tile = provide_STATVAR(tile)
            
        end
        
        
        %assemble the stratigraphy
        function tile = finalize_init(tile)
            
        end
        
        
        function tile = interpolate_forcing_tile(tile)
            
        end
        
        
        
        function tile = store_OUT_tile(tile)
            
        end
        
        
        
        function tile = run_model(tile)
            disp('rearranging forcing files step 2')
            
            number_of_cores = tile.RUN_INFO.PARA.num_ranks;
            my_core = tile.RUN_INFO.PARA.my_rank;
            
            file_folder = [tile.PARA.file_path tile.RUN_INFO.PARA.run_name '/'];
            run_name =tile.RUN_INFO.PARA.run_name;
            variables = {'ERA_melt_bare'; 'ERA_melt_forest'; 'ERA_snowfall_downscaled'; 'ERA_T_downscaled'; 'final_av_T'; 'final_MODIS_weight'; 'timestamp'};
            
            for i=1:size(variables,1)-1
                forcing.DATA.(variables{i,1}) = [];
                for j=1:number_of_cores
                    data = ncread([file_folder tile.PARA.in_name num2str(j) '_' num2str(my_core) '.nc'], variables{i,1}, [1 1], [Inf Inf], [1 1]);
                    forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}) data];
                end
                
                
                %change filename
%                 nccreate([file_folder run_name  '_' num2str(my_core) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
%                     {'x',size(forcing.DATA.(variables{i,1}),1),'y',size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);
%                 
%                 ncwrite([file_folder run_name  '_' num2str(my_core) '.nc'], ...
%                     variables{i,1}, forcing.DATA.(variables{i,1}), [1 1], [1 1]);

                nccreate([file_folder run_name '_' variables{i,1} '_' num2str(my_core) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
                    {'x',size(forcing.DATA.(variables{i,1}),1),'y',size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);
                
                ncwrite([file_folder run_name '_' variables{i,1} '_' num2str(my_core) '.nc'], ...
                    variables{i,1}, forcing.DATA.(variables{i,1}), [1 1], [1 1]);
                
                
            end
            
            %timestamp
            i=size(variables,1);
            forcing.DATA.(variables{i,1}) = [];
            for j=1:number_of_cores
                data = ncread([file_folder tile.PARA.in_name num2str(j) '_' num2str(my_core) '.nc'], variables{i,1}, [1], [Inf], [ 1]);
                forcing.DATA.(variables{i,1}) = [forcing.DATA.(variables{i,1}) data'];
            end
            
            
            %change filename
            nccreate([file_folder run_name '_' variables{i,1} '_' num2str(my_core) '.nc'], variables{i,1}, 'datatype', 'single', 'Dimensions', ...
                {'y',size(forcing.DATA.(variables{i,1}),2)}, 'FillValue', -9999);
            
            ncwrite([file_folder run_name '_' variables{i,1} '_' num2str(my_core) '.nc'], ...
                variables{i,1}, forcing.DATA.(variables{i,1}), [1], [1]);

            
        end
        

    end
end



