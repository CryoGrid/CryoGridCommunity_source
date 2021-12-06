% base class build a model tile

classdef TILE_3D_singleClass < matlab.mixin.Copyable
    
    properties
        

        BUILDER
        PARA
        RUN_INFO
        FORCING
        CONST
        GRID
        OUT        
        
        t        
        timestep
        next_break_time
                
        LATERAL
        GROUND

    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)

            tile.PARA.builder = [];
            
            %new_init
            tile.PARA.latitude = [];
            tile.PARA.longitude = [];
            tile.PARA.altitude = [];
            tile.PARA.domain_depth = [];
            tile.PARA.area = [];
            
            tile.PARA.forcing_class = [];
            tile.PARA.forcing_class_index = [];
            tile.PARA.grid_class = [];
            tile.PARA.grid_class_index = [];
            tile.PARA.out_class = [];
            tile.PARA.out_class_index = [];
            tile.PARA.ground_class = [];  %single class
            tile.PARA.ground_class_index = [];
            tile.PARA.strat_statvar_class = [];
            tile.PARA.strat_statvar_class_index = [];
%             tile.PARA.lateral_class = [];
%             tile.PARA.lateral_class_index = [];
%             tile.PARA.lateral_IA_classes = [];
%             tile.PARA.lateral_IA_classes_index = [];
            
            %restart_OUT_last_timestep
            tile.PARA.restart_file_path = [];
            tile.PARA.restart_file_name = [];
            
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

        function tile = interact_lateral(tile)
            tile.LATERAL = interact(tile.LATERAL, tile);
        end
        
        function tile = store_OUT_tile(tile)
            tile.OUT = store_OUT(tile.OUT, tile);
        end        
        
        
        
        function tile = run_model(tile)
            
%             TOP_CLASS = tile.TOP_CLASS;
%             BOTTOM_CLASS = tile.BOTTOM_CLASS;
            GROUND = tile.GROUND;

            %=========================================================================
            %TIME INTEGRATION
            %=========================================================================
            while tile.t < tile.FORCING.PARA.end_time
                
                
                %interpolate focing data to time t
                tile = interpolate_forcing_tile(tile);
                
                %upper boundar condition (uppermost class only)
                GROUND = get_boundary_condition_u(GROUND, tile);
                
                
                %lower boundary condition (lowermost class)
                GROUND = get_boundary_condition_l(GROUND,  tile);  %At this point, CURRENT is equal to BOTTOM_CLASS
                
                %calculate spatial derivatives
                GROUND = get_derivatives_prognostic(GROUND, tile);

                %calculate timestep [second]
                tile.timestep = min(tile.timestep, get_timestep(GROUND, tile));
                
                %prognostic step - integrate prognostic variables in time
                GROUND = advance_prognostic(GROUND, tile);
                
                %diagnostic step - compute diagnostic variables
                GROUND = compute_diagnostic_first_cell(GROUND, tile); %calculate Lstar, only uppermost class

                GROUND = compute_diagnostic(GROUND, tile);

                %triggers
                GROUND = check_trigger(GROUND, tile);

                %tile = interact_lateral(tile);

                %update time variable t
                tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
                
                %model
                tile = store_OUT_tile(tile);
            end
            
        end
        
        
        %---BUILDER functions--------------

        function check_if_PARA_assigned(tile)
            if size(tile.PARA.builder,1) == 0 && size(tile.PARA.builder,2) == 0
                disp(['PARA builder in class ' class(tile) ' not assigned'])
            end
            if strcmp(tile.PARA.builder, 'new_init')
                parameters = {'latitude'; 'longitude'; 'altitude'; 'domain_depth'; 'area'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
                            'out_class_index'; 'ground_class'; 'ground_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'};
            elseif strcmp(tile.PARA.builder, 'update_forcing_out')
                parameters = { 'forcing_class'; 'forcing_class_index';  'out_class'; 'out_class_index'};
            elseif strcmp(tile.PARA.builder, 'update_forcing_out')
                parameters = {'restart_file_path'; 'restart_file_name'};
            else
                parameters = fieldnames(tile.PARA);
            end
            for i=1:size(parameters,1)
                if size(tile.PARA.(parameters{i,1}),1) == 0 && size(tile.PARA.(parameters{i,1}),2) == 0
                    disp(['Warning: PARA ' parameters{i,1} ' in class ' class(tile) ' not assigned'])
                end
            end
        end
        
        function tile = build_tile_new_init(tile)
            
            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
            
            %1. forcing
            %tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.FORCING{tile.PARA.forcing_index,1});
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
            tile.FORCING = finalize_init(tile.FORCING, tile);
            
            %2. grid
            %tile.GRID = copy(tile.RUN_INFO.PPROVIDER.FUNCTIONAL_CLASSES.GRID{tile.PARA.grid_index,1});
            tile.GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.grid_class){tile.PARA.grid_class_index,1});
            tile.GRID = finalize_init(tile.GRID, tile);
            
            %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes 
            for i=1:size(tile.PARA.strat_statvar_class,1)
                strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_statvar_class{i,1}){tile.PARA.strat_statvar_class_index(i,1),1});
                strat_statvar_class = finalize_init(strat_statvar_class, tile);
            end

            %4. build stratigraphy
            tile.GROUND = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.ground_class){tile.PARA.ground_class,1};

            %5. assign STATVAR using STRATGRAPHY_STATVAR classes            
            variables = fieldnames(GROUND.STATVAR);
            for j=1:size(variables,1)
                if isfield(tile.GRID.STATVAR, variables{j,1})
                    GROUND.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1});
                end
            end

            %6. finalize initialization for subsurface classes
            GROUND = finalize_init(GROUND, tile);
            
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            
%             %11. assign LATERAL classes 
%             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
%             tile.LATERAL = finalize_init(tile.LATERAL, tile);
            
            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
           
        end
        
%         function tile = build_tile_update_forcing_out(tile)
% 
%             
%             %2. forcing -> special forcing class required
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
%             
%             %use old tile
%             tile.GRID = tile.RUN_INFO.TILE.GRID;   
%             tile.BOTTOM = tile.RUN_INFO.TILE.BOTTOM;
%             tile.TOP = tile.RUN_INFO.TILE.TOP;
%             tile.BOTTOM_CLASS = tile.RUN_INFO.TILE.BOTTOM_CLASS;
%             tile.TOP_CLASS = tile.RUN_INFO.TILE.TOP_CLASS;
%             tile.timestep = tile.RUN_INFO.TILE.timestep;
%             tile.LATERAL = tile.RUN_INFO.TILE.LATERAL;            
%             
%             %use old PARA, but overwrite all newly set values
%             PARA_new = tile.PARA;
%             tile.PARA = tile.RUN_INFO.TILE.PARA;
%             fn = fieldnames(PARA_new);
%             for i=1:size(fn,2)
%                 tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
%             end
% 
%             
%             tile.FORCING = finalize_init(tile.FORCING, tile); 
%             tile.OUT = finalize_init(tile.OUT, tile);           
%             %10. assign time, etc.
%             tile.t = tile.FORCING.PARA.start_time;
%             
%             %reset IA time
%             tile.LATERAL.IA_TIME = tile.t + tile.LATERAL.IA_TIME_INCREMENT;
%             
%             %reset time for BGC class (do mothing if no BGC class exists)
%             CURRENT = tile.TOP.NEXT;
%             while ~isequal(CURRENT.NEXT, tile.BOTTOM)
%                 CURRENT = reset_time_BGC(CURRENT, tile);
%                 CURRENT = CURRENT.NEXT;
%             end
%             
% %             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
% %             tile.LATERAL = finalize_init(tile.LATERAL, tile);
%             
%             tile.RUN_INFO.TILE = tile;
% 
%             
%             tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
%             tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
%         end
%         
%         
%         function tile = build_tile_restart_OUT_last_timestep(tile)
%             temp=load([tile.PARA.restart_file_path tile.PARA.restart_file_name]);
%             variables = fieldnames(temp.out.STRATIGRAPHY);
%             for i=1:size(variables,1)
%                 tile.(variables{i,1}) = temp.out.STRATIGRAPHY.(variables{i,1});
%             end
% %             tile.LATERAL.IA_CLASSES = {};
% %             tile.LATERAL.PARA.num_realizations = 1;
% %             tile.LATERAL.PARA.worker_number = 1;
% %             tile.OUT.OUTPUT_TIME = tile.OUT.OUTPUT_TIME+100;
%             %tile.LATERAL.IA_TIME = tile.FORCING.PARA.end_time;
%         end
%         
        

    end
end



