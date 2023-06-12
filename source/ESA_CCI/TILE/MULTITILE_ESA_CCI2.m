
classdef MULTITILE_ESA_CCI2 < matlab.mixin.Copyable
    
    properties
        PARA
        RUN_INFO
        BUILDER
        FORCING
        CONST
        GRID
        OUT        
        SUBSURFACE_CLASS
        ENSEMBLE
        STRATIGRAPHY
        
        t        
        timestep
        
        TOP
        BOTTOM
    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)
            
            tile.PARA.builder = [];
            tile.PARA.domain_depth = [];
            tile.PARA.number_of_realizations = [];
            tile.PARA.ensemble_size = [];
            tile.PARA.ensemble_class =[];
            tile.PARA.ensemble_class_index =[];
            tile.PARA.subsurface_class =[];
            tile.PARA.subsurface_class_index =[];
            tile.PARA.forcing_class = [];
            tile.PARA.forcing_class_index = [];
            tile.PARA.strat_statvar_class = [];
            tile.PARA.strat_statvar_class_index = [];
%             tile.PARA.strat_statvar_class_index_end = [];

            tile.PARA.grid_class = [];
            tile.PARA.grid_class_index = [];
            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];
            
        end
        
        function tile = provide_CONST(tile)

            tile.CONST.day_sec = [];

        end
        
        function tile = provide_STATVAR(tile)

        end
        

       
        
        %assemble the stratigraphy
        function tile = finalize_init(tile)
            
            builder = str2func(tile.PARA.builder);
            tile.BUILDER = builder();
            tile.BUILDER.TILE = tile;
            
            build_tile(tile.BUILDER);
           
        end
        

        function tile = interpolate_forcing_tile(tile)
             tile.FORCING = interpolate_forcing(tile.FORCING, tile);
        end

        
        function tile = store_OUT_tile(tile)
            for i=1:size(tile.OUT, 1)
                tile.OUT{i,1} = store_OUT(tile.OUT{i,1}, tile);
            end
        end  
        
        
        
        function tile = run_model(tile)
            
            %=========================================================================
            %TIME INTEGRATION
            %=========================================================================
            while tile.t < tile.FORCING.PARA.end_time
                
                CURRENT = tile.SUBSURFACE_CLASS;
                
                %interpolate focing data to time t
                tile = interpolate_forcing_tile(tile);
                
                %upper boundar condition (uppermost class only)
                CURRENT = get_boundary_condition_u(CURRENT, tile);
                
                %lower boundary condition (lowermost class)
                CURRENT = get_boundary_condition_l(CURRENT,  tile);
                
                %calculate spatial derivatives
                CURRENT = get_derivatives_prognostic(CURRENT, tile);

                %prognostic step - integrate prognostic variables in time
                CURRENT = advance_prognostic(CURRENT, tile);
                
                %diagnostic step - compute diagnostic variables
                CURRENT = compute_diagnostic(CURRENT, tile);

                %triggers
                CURRENT = check_trigger(CURRENT, tile);

                %update time variable t
                tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
                
                %model output
                tile = store_OUT_tile(tile);
            end
            
        end
        
        
        %---BUILDER functions--------------
        
        function tile = build_tile_new_init(tile)
            
            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            tile.PARA.tile_size = tile.PARA.number_of_realizations .* tile.PARA.ensemble_size;

            %assign ensemble classes
%             for i = 1:size(tile.PARA.ensemble_class,1)
%               tile.ENSEMBLE{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class{i,1}){tile.PARA.ensemble_class_index(i,1)});  
%               tile.ENSEMBLE{i,1} = finalize_init(tile.ENSEMBLE{i,1}, tile);  
%             end
            tile.ENSEMBLE = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class){tile.PARA.ensemble_class_index});
            tile.ENSEMBLE = finalize_init(tile.ENSEMBLE, tile);
                        
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);


            %3. grid -> existing GRID class should do, possibly some
            %add-ons for snow required
            tile.GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.grid_class){tile.PARA.grid_class_index,1});
            tile.GRID = finalize_init(tile.GRID, tile);
                        
            %done by GRID class
%             tile.GRID.STATVAR.depth = repmat(tile.GRID.STATVAR.GRID, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.midPoints = repmat(tile.GRID.STATVAR.MIDPOINTS, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.layerThick = repmat(tile.GRID.STATVAR.layerThick, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.layerDistance = repmat(tile.GRID.STATVAR.layerDistance, 1, tile.PARA.tile_size);

            %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
            for i=1:size(tile.PARA.strat_statvar_class,1)
                strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
                strat_statvar_class = finalize_init(strat_statvar_class, tile);
            end

%             %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
%             for i=1:size(tile.PARA.strat_statvar_class,1) %this is only a single class
%                 for j = tile.PARA.strat_statvar_class_index_start(i,1):tile.PARA.strat_statvar_class_index_end(i,1)
%                     tile.STRATIGRAPHY{i,j} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){j,1});
%                     tile.STRATIGRAPHY{i,j} = finalize_init(tile.STRATIGRAPHY{i,j}, tile); 
%                 end
%             end
            
%             tile.ENSEMBLE = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class){tile.PARA.ensemble_class_index,1});
%             tile.ENSEMBLE = finalize_init(tile.ENSEMBLE, tile);
%             ensemble = tile.ENSEMBLE.STATVAR;
%             save([tile.RUN_INFO.PPROVIDER.PARA.result_path tile.RUN_INFO.PPROVIDER.PARA.run_name '/' tile.RUN_INFO.PARA.run_name '_ensemble.mat'], 'ensemble');
%           
            %1. build stratigraphy
            tile.SUBSURFACE_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.subsurface_class){tile.PARA.subsurface_class_index,1});
            variables = fieldnames(tile.SUBSURFACE_CLASS.STATVAR);
            for j=1:size(variables,1)
                if isfield(tile.GRID.STATVAR, variables{j,1})
                    tile.SUBSURFACE_CLASS.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1});
                end
            end
       
            tile.SUBSURFACE_CLASS = finalize_init(tile.SUBSURFACE_CLASS, tile);
            
            %set this for being able to use existing OUT classes
            tile.TOP = Top();
            tile.BOTTOM = Bottom();
            tile.TOP.NEXT = tile.SUBSURFACE_CLASS;
            tile.TOP.NEXT.PREVIOUS = tile.TOP;
            tile.SUBSURFACE_CLASS.NEXT = tile.BOTTOM();
            tile.BOTTOM.PREVIOUS = tile.SUBSURFACE_CLASS;            
            
            %set fixed timestep!
            tile.timestep = tile.SUBSURFACE_CLASS.PARA.timestep;
%             
%             %2. forcing -> special forcing class required
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
%             tile.FORCING = finalize_init(tile.FORCING, tile);
% 
             %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            %this is not solved in a good way
            tile.SUBSURFACE_CLASS.TEMP.adjust_stratigraphy_date = datenum([tile.SUBSURFACE_CLASS.PARA.adjust_stratigraphy_date datestr(tile.t, 'yyyy')], 'dd.mm.yyyy');

            %12. assign OUT classes
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class{i,1}){tile.PARA.out_class_index(i,1),1});
                tile.OUT{i,1} = finalize_init(tile.OUT{i,1}, tile);
            end
            
            tile.RUN_INFO.TILE = tile;
            
        end
        
        

        function tile = build_tile_new_init_TTOP_GlobPermafrost(tile)
            
                       
            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            %this must be done by assign Tile properties class
            tile.PARA.tile_size = tile.PARA.number_of_realizations .* tile.PARA.ensemble_size;

            %assign ensemble classes
%             for i = 1:size(tile.PARA.ensemble_class,1)
%               tile.ENSEMBLE{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class{i,1}){tile.PARA.ensemble_class_index(i,1)});  
%               tile.ENSEMBLE{i,1} = finalize_init(tile.ENSEMBLE{i,1}, tile);  
%             end
            tile.ENSEMBLE = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ensemble_class){tile.PARA.ensemble_class_index});
            tile.ENSEMBLE = finalize_init(tile.ENSEMBLE, tile);
            
            %3. grid -> existing GRID class should do, possibly some
            %add-ons for snow required
            tile.GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.grid_class){tile.PARA.grid_class_index,1});
            tile.GRID = finalize_init(tile.GRID, tile);
                        
            %done by GRID class
%             tile.GRID.STATVAR.depth = repmat(tile.GRID.STATVAR.GRID, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.midPoints = repmat(tile.GRID.STATVAR.MIDPOINTS, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.layerThick = repmat(tile.GRID.STATVAR.layerThick, 1, tile.PARA.tile_size);
%             tile.GRID.STATVAR.layerDistance = repmat(tile.GRID.STATVAR.layerDistance, 1, tile.PARA.tile_size);

            %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
            for i=1:size(tile.PARA.strat_statvar_class,1)
                strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
                strat_statvar_class = finalize_init(strat_statvar_class, tile);
            end

            tile.timestep = 1; %random value, so that data can be loaded
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);
            
            TTOP_Globpermafrost = copy(tile.RUN_INFO.PPROVIDER.CLASSES.TTOP_MULTITILE_GlobPF{1,1}); %must be made a parameter!
            tile.GRID.STATVAR.T = get_TTOP_from_forcing(TTOP_Globpermafrost, tile);      
            tile.PARA.TTOP = TTOP_Globpermafrost;
            
            %1. build stratigraphy
            tile.SUBSURFACE_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.subsurface_class){tile.PARA.subsurface_class_index,1});
            variables = fieldnames(tile.SUBSURFACE_CLASS.STATVAR);
            for j=1:size(variables,1)
                if isfield(tile.GRID.STATVAR, variables{j,1})
                    tile.SUBSURFACE_CLASS.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1});
                end
            end

            tile.SUBSURFACE_CLASS = finalize_init(tile.SUBSURFACE_CLASS, tile);
            
            %set this for being able to use existing OUT classes
            tile.TOP = Top();
            tile.BOTTOM = Bottom();
            tile.TOP.NEXT = tile.SUBSURFACE_CLASS;
            tile.TOP.NEXT.PREVIOUS = tile.TOP;
            tile.SUBSURFACE_CLASS.NEXT = tile.BOTTOM();
            tile.BOTTOM.PREVIOUS = tile.SUBSURFACE_CLASS;            
            

            %set fixed timestep!
            tile.timestep = tile.SUBSURFACE_CLASS.PARA.timestep;
            
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING.PARA.preprocessed = 1;
            tile.FORCING = finalize_init(tile.FORCING, tile);

             %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            %this is not solved in a good way
            tile.SUBSURFACE_CLASS.TEMP.adjust_stratigraphy_date = datenum([tile.SUBSURFACE_CLASS.PARA.adjust_stratigraphy_date datestr(tile.t, 'yyyy')], 'dd.mm.yyyy');

            %12. assign OUT classes
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class{i,1}){tile.PARA.out_class_index(i,1),1});
                tile.OUT{i,1} = finalize_init(tile.OUT{i,1}, tile);
            end
            
            tile.RUN_INFO.TILE = tile;
            
        end
        

        function tile = build_tile_update_forcing_out(tile)
                        
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});

            %12. assign OUT classes
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class{i,1}){tile.PARA.out_class_index(i,1),1});
            end
            
            %use old tile
            tile.GRID = tile.RUN_INFO.TILE.GRID;   
            tile.ENSEMBLE = tile.RUN_INFO.TILE.ENSEMBLE;
            tile.SUBSURFACE_CLASS = tile.RUN_INFO.TILE.SUBSURFACE_CLASS;
            tile.BOTTOM = tile.RUN_INFO.TILE.BOTTOM;
            tile.TOP = tile.RUN_INFO.TILE.TOP;
            tile.timestep = tile.RUN_INFO.TILE.timestep;
            
            %use old PARA, but overwrite all newly set values
            PARA_new = tile.PARA;
            tile.PARA = tile.RUN_INFO.TILE.PARA;
            fn = fieldnames(PARA_new);
            for i=1:size(fn,1)
                if ~isempty(PARA_new.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
                end
            end

            
            tile.FORCING = finalize_init(tile.FORCING, tile); 
            
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = finalize_init(tile.OUT{i,1}, tile);
            end
            
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            
            tile.RUN_INFO.TILE = tile;
            
        end
        
        function tile = build_tile_update_forcing_out_TTOP(tile)

            rk = mean(tile.RUN_INFO.TILE.SUBSURFACE_CLASS.STATVAR.k_thawed(1:10,:),1)./ mean(tile.RUN_INFO.TILE.SUBSURFACE_CLASS.STATVAR.k_frozen(1:10,:),1); %uppermost meter as proxy for AL/seasonal frost layer
            TTOP_Globpermafrost = copy(tile.RUN_INFO.PPROVIDER.CLASSES.TTOP_MULTITILE_GlobPF{1,1});
            TTOP = get_TTOP(TTOP_Globpermafrost, tile.RUN_INFO.TILE.OUT.TEMP.FDD, tile.RUN_INFO.TILE.OUT.TEMP.TDD, rk, tile.RUN_INFO.TILE.OUT.TEMP.N);
      
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            
            %12. assign OUT classes
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class{i,1}){tile.PARA.out_class_index(i,1),1});
            end
            
            %use old tile
            tile.GRID = tile.RUN_INFO.TILE.GRID;
            tile.GRID.STATVAR.T = TTOP; %assign new TTOP
            tile.ENSEMBLE = tile.RUN_INFO.TILE.ENSEMBLE;
            tile.timestep = tile.RUN_INFO.TILE.timestep;
            
            %use old PARA, but overwrite all newly set values
            PARA_new = tile.PARA;
            tile.PARA = tile.RUN_INFO.TILE.PARA;
            fn = fieldnames(PARA_new);
            for i=1:size(fn,1)
                if ~isempty(PARA_new.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
                end
            end
            
            
            %1. build stratigraphy
            tile.SUBSURFACE_CLASS = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.subsurface_class){tile.PARA.subsurface_class_index,1});
            variables = fieldnames(tile.SUBSURFACE_CLASS.STATVAR);
            for j=1:size(variables,1)
                if isfield(tile.GRID.STATVAR, variables{j,1})
                    tile.SUBSURFACE_CLASS.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1});
                end
            end
            %use water contents from previous run 
            tile.SUBSURFACE_CLASS.STATVAR.waterIce = tile.RUN_INFO.TILE.SUBSURFACE_CLASS.STATVAR.waterIce ./ tile.RUN_INFO.TILE.SUBSURFACE_CLASS.STATVAR.layerThick(5:end,:);

            tile.SUBSURFACE_CLASS = finalize_init(tile.SUBSURFACE_CLASS, tile);
            
            %set this for being able to use existing OUT classes
            tile.TOP = Top();
            tile.BOTTOM = Bottom();
            tile.TOP.NEXT = tile.SUBSURFACE_CLASS;
            tile.TOP.NEXT.PREVIOUS = tile.TOP;
            tile.SUBSURFACE_CLASS.NEXT = tile.BOTTOM();
            tile.BOTTOM.PREVIOUS = tile.SUBSURFACE_CLASS;            
            

            %set fixed timestep!
            tile.timestep = tile.SUBSURFACE_CLASS.PARA.timestep;
            
            %2. forcing -> special forcing class required
            tile.FORCING = finalize_init(tile.FORCING, tile);

             %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            %this is not solved in a good way
       %     tile.SUBSURFACE_CLASS.TEMP.adjust_stratigraphy_date = datenum([tile.SUBSURFACE_CLASS.PARA.adjust_stratigraphy_date datestr(tile.t, 'yyyy')], 'dd.mm.yyyy');

            %12. assign OUT classes
            for i=1:size(tile.PARA.out_class,1)
                tile.OUT{i,1} = finalize_init(tile.OUT{i,1}, tile);
            end
            
            tile.RUN_INFO.TILE = tile;
            
        end
        
        

    end
end



