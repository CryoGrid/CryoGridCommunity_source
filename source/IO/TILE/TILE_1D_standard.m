%========================================================================
% CryoGrid TILE class TILE_1D_standard
% TILE class designed for multi-physics simulations with a stratigraphy of
% connected CryoGrid stratigraphy classes. TILE_1D_standard represents all 
% aspects of a single model simulation, with model forcing, storage of
% model output, etc. It supports both lateral coupling to external
% environments and between different TILE_1D_standard classes in a parallel
% fashion
% The initialization procedure is determined by the choice of the 
% TILE_BUILDER class (set by PARA "builder") which also determines the 
% paramter set which must be defined by the users. 

% S. Westermann, Dec 2020
%========================================================================

classdef TILE_1D_standard < matlab.mixin.Copyable
    
    properties
        
        BUILDER
        PARA
        RUN_INFO
        FORCING
        CONST
        GRID
        OUT        
        STORE
        TEMP
        
        t        
        timestep
        next_break_time
                
        LATERAL
        TOP          % Indicates current top position
        BOTTOM       % Indicates current bottom position
        TOP_CLASS    % Uppermost ground/snow class in the soil column
        BOTTOM_CLASS % Lowermost ground class in the soil column

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
            tile.PARA.strat_classes_class = [];
            tile.PARA.strat_classes_class_index = [];
            tile.PARA.strat_statvar_class = [];
            tile.PARA.strat_statvar_class_index = [];
            tile.PARA.lateral_class = [];
            tile.PARA.lateral_class_index = [];
            tile.PARA.lateral_IA_classes = [];
            tile.PARA.lateral_IA_classes_index = [];
            
            %new_init_steady_state
            tile.PARA.init_steady_state_class = []; % if left empty, use T_first_cell
            tile.PARA.init_steady_state_class_index = [];
            tile.PARA.T_first_cell = [];
            tile.PARA.start_depth_steady_state = [];
            
            %restart_OUT_last_timestep
            tile.PARA.restart_file_path = [];
            tile.PARA.restart_file_name = [];

            tile.PARA.unit_conversion_class = 'UNIT_CONVERSION_standard'; %can be overwritten if needed
            
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
            TOP = tile.TOP;
            BOTTOM = tile.BOTTOM;
            TOP.LATERAL = tile.LATERAL;

            %=========================================================================
            %TIME INTEGRATION
            %=========================================================================
            while tile.t < tile.FORCING.PARA.end_time
                
                %interpolate focing data to time t
                tile = interpolate_forcing_tile(tile);
                
                %upper boundar condition (uppermost class only)
                TOP.NEXT = get_boundary_condition_u(TOP.NEXT, tile);
                
                %set fluxes between classes in the stratigrapht
                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT.NEXT, BOTTOM)
                    get_boundary_condition_m(CURRENT.IA_NEXT, tile); %call interaction class function
                    CURRENT = CURRENT.NEXT;
                end
                
                %lower boundary condition (lowermost class)
                CURRENT = get_boundary_condition_l(CURRENT,  tile);  %At this point, CURRENT is equal to BOTTOM_CLASS
                
                %calculate spatial derivatives
                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT, BOTTOM)
                    CURRENT = get_derivatives_prognostic(CURRENT, tile);
                    CURRENT = CURRENT.NEXT;
                end
                
                %calculate timestep [second]
                CURRENT = TOP.NEXT;
                tile.timestep = 1e8;
                while ~isequal(CURRENT, BOTTOM)
                    tile.timestep = min(tile.timestep, get_timestep(CURRENT, tile));
                    CURRENT = CURRENT.NEXT;
                end
                tile.next_break_time = min(tile.LATERAL.IA_TIME, tile.OUT.OUTPUT_TIME);
                tile.timestep = min(tile.timestep, (tile.next_break_time - tile.t).*tile.CONST.day_sec);
                
                %prognostic step - integrate prognostic variables in time
                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT, BOTTOM)
                    CURRENT = advance_prognostic(CURRENT, tile);
                    CURRENT = CURRENT.NEXT;
                end
                
                %diagnostic step - compute diagnostic variables
                TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, tile); %calculate Lstar, only uppermost class
                CURRENT = BOTTOM.PREVIOUS;
                while ~isequal(CURRENT, TOP)
                    CURRENT = compute_diagnostic(CURRENT, tile);
                    CURRENT = CURRENT.PREVIOUS;
                end
                
                
                %triggers
                CURRENT = TOP.NEXT;
                while ~isequal(CURRENT, BOTTOM)
                    CURRENT = check_trigger(CURRENT, tile);
                    CURRENT = CURRENT.NEXT;
                end

                %lateral interactions
                tile = interact_lateral(tile);
                                
                %set TOP_CLASS and BOTTOM_CLASS for convenient access
                tile.TOP_CLASS = TOP.NEXT;
                tile.BOTTOM_CLASS = BOTTOM.PREVIOUS;
                
                %update time variable t
                tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
                
                %model
                tile = store_OUT_tile(tile);
            end

            % extra carriage return needed for OUT classes that use
            % fprintf to print and overwrite step date output to 
            % console window.
            fprintf('\n\n')
            
        end
        
        
        %---BUILDER functions--------------

        function check_if_PARA_assigned(tile)
            if size(tile.PARA.builder,1) == 0 && size(tile.PARA.builder,2) == 0
                disp(['PARA builder in class ' class(tile) ' not assigned'])
            end
            if strcmp(tile.PARA.builder, 'new_init') 
                parameters = {'latitude'; 'longitude'; 'altitude'; 'domain_depth'; 'area'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
                            'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
                            'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'};
            elseif strcmp(tile.PARA.builder, 'new_init_steady_state')
                parameters = {'latitude'; 'longitude'; 'altitude'; 'domain_depth'; 'area'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
                    'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
                    'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'; 'init_steady_state_class'; 'init_steady_state_class_index';...
                    'T_first_cell'; 'start_depth_steady_state'};
            elseif strcmp(tile.PARA.builder, 'update_forcing_out')
                parameters = { 'forcing_class'; 'forcing_class_index';  'out_class'; 'out_class_index'};
            elseif strcmp(tile.PARA.builder, 'restart_OUT_last_timestep')
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
            class_list = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.class_name;
            class_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.class_index;
            tile.TOP = Top();
            CURRENT = tile.TOP;
            for i=1:size(class_list,1)
                CURRENT.NEXT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(class_list{i,1}){class_index(i,1)});
                CURRENT.NEXT.PREVIOUS = CURRENT;
                CURRENT = CURRENT.NEXT;
            end
            tile.BOTTOM = Bottom();
            CURRENT.NEXT = tile.BOTTOM;
            tile.BOTTOM.PREVIOUS = CURRENT;
            
            tile.TOP_CLASS = tile.TOP.NEXT;
            tile.BOTTOM_CLASS = tile.BOTTOM.PREVIOUS;
            
            %5. assign STATVAR using STRATGRAPHY_STATVAR classes
            class_depths = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.depth;
            class_depths = [class_depths; tile.GRID.STATVAR.GRID(end,1)];
            
            CURRENT = tile.TOP_CLASS;
            for i=1:size(class_list,1)
                variables = fieldnames(CURRENT.STATVAR);
                range = (tile.GRID.STATVAR.MIDPOINTS > class_depths(i,1) & tile.GRID.STATVAR.MIDPOINTS <= class_depths(i+1,1));
                %CURRENT.STATVAR.layerThick = tile.GRID.STATVAR.LAYERTHICK(range,1);
                CURRENT.STATVAR.upperPos = tile.PARA.altitude - class_depths(i,1);
                CURRENT.STATVAR.lowerPos = tile.PARA.altitude - class_depths(i+1,1);
                for j=1:size(variables,1)
                    if isfield(tile.GRID.STATVAR, variables{j,1})
                        CURRENT.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1})(range);
                    end
                end
                CURRENT=CURRENT.NEXT;
            end
            
            %6. set top depths relative to surface and finalize initialization for
            %subsurface classes
            
            CURRENT = tile.TOP_CLASS;
            CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
            
            CURRENT = convert_units(CURRENT, tile);
            CURRENT = finalize_init(CURRENT, tile);

            CURRENT.PARA.target_grid = tile.GRID.STATVAR.GRID;
            CURRENT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
            while ~isequal(CURRENT.NEXT, tile.BOTTOM_CLASS.NEXT)
                CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
                
                CURRENT.NEXT = convert_units(CURRENT.NEXT, tile);
                CURRENT.NEXT = finalize_init(CURRENT.NEXT, tile);

                CURRENT.NEXT.PARA.target_grid = tile.GRID.STATVAR.GRID;
                CURRENT.NEXT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
                
                CURRENT = CURRENT.NEXT;
            end

            %7. assign interaction classes 
            CURRENT = tile.TOP_CLASS;
            
            while ~isequal(CURRENT.NEXT, tile.BOTTOM)
                ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
                CURRENT.IA_NEXT = ia_class;
                CURRENT.IA_NEXT.PREVIOUS = CURRENT;
                CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
                CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
                
                finalize_init(CURRENT.IA_NEXT, tile); %Added Sebastian, defined in IA_BASE as empty
                
                CURRENT = CURRENT.NEXT;
            end
            
            %8. assign SNOW class
            snow_class_name = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_name;
            snow_class_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_index;
            
            if ~isempty(snow_class_name) && sum(isnan(snow_class_name))==0
                snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
                snow_class = snow_class{snow_class_index,1};
                
                %tile.TOP.STORE.SNOW = copy(snow_class);
                tile.STORE.SNOW = copy(snow_class);
                %no convert_units for SNOW classes needed at this point!
                %tile.TOP.STORE.SNOW = finalize_init(tile.TOP.STORE.SNOW, tile); %make this dependent on TILE!
                tile.STORE.SNOW = finalize_init(tile.STORE.SNOW, tile); %make this dependent on TILE!
            end
            
            %9. assign sleeping classes
            sleeping_classes = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_name;
            sleeping_classes_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_index; 

%             for i=1:size(sleeping_classes,1)
%                 sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
%                 sc = sc{sleeping_classes_index(i,1),1};
%                 tile.TOP.STORE.SLEEPING{i,1} = copy(sc);
%                 tile.TOP.STORE.SLEEPING{i,1} = convert_units(tile.TOP.STORE.SLEEPING{i,1}, tile);
%                 tile.TOP.STORE.SLEEPING{i,1} = finalize_init(tile.TOP.STORE.SLEEPING{i,1}, tile);
%                 tile.TOP.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
%             end
            for i=1:size(sleeping_classes,1)
                sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
                sc = sc{sleeping_classes_index(i,1),1};
                tile.STORE.SLEEPING{i,1} = copy(sc);
                tile.STORE.SLEEPING{i,1} = convert_units(tile.STORE.SLEEPING{i,1}, tile);
                tile.STORE.SLEEPING{i,1} = finalize_init(tile.STORE.SLEEPING{i,1}, tile);
                tile.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
            end
            
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            
            %11. assign LATERAL classes 
            tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
            tile.LATERAL = finalize_init(tile.LATERAL, tile);
            
            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
           
        end
        
        function tile = build_tile_new_init_steady_state(tile)
            
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
            class_list = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.class_name;
            class_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.class_index;
            tile.TOP = Top();
            CURRENT = tile.TOP;
            for i=1:size(class_list,1)
                CURRENT.NEXT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(class_list{i,1}){class_index(i,1)});
                CURRENT.NEXT.PREVIOUS = CURRENT;
                CURRENT = CURRENT.NEXT;
            end
            tile.BOTTOM = Bottom();
            CURRENT.NEXT = tile.BOTTOM;
            tile.BOTTOM.PREVIOUS = CURRENT;
            
            tile.TOP_CLASS = tile.TOP.NEXT;
            tile.BOTTOM_CLASS = tile.BOTTOM.PREVIOUS;
            
            %5. assign STATVAR using STRATGRAPHY_STATVAR classes
            class_depths = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.classes.depth;
            class_depths = [class_depths; tile.GRID.STATVAR.GRID(end,1)];
            
            CURRENT = tile.TOP_CLASS;
            for i=1:size(class_list,1)
                variables = fieldnames(CURRENT.STATVAR);
                range = (tile.GRID.STATVAR.MIDPOINTS > class_depths(i,1) & tile.GRID.STATVAR.MIDPOINTS <= class_depths(i+1,1));
                %CURRENT.STATVAR.layerThick = tile.GRID.STATVAR.LAYERTHICK(range,1);
                CURRENT.STATVAR.upperPos = tile.PARA.altitude - class_depths(i,1);
                CURRENT.STATVAR.lowerPos = tile.PARA.altitude - class_depths(i+1,1);
                for j=1:size(variables,1)
                    if isfield(tile.GRID.STATVAR, variables{j,1})
                        CURRENT.STATVAR.(variables{j,1}) = tile.GRID.STATVAR.(variables{j,1})(range);
                    end
                end
                CURRENT=CURRENT.NEXT;
            end
            
            %6. set top depths relative to surface and finalize initialization for
            %subsurface classes
            
            %determine T_first_cell from class 
            if ~isempty(tile.PARA.init_steady_state_class) && sum(isnan(tile.PARA.init_steady_state_class)) == 0
                init_steady_state_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.init_steady_state_class){tile.PARA.init_steady_state_class_index,1});
                init_steady_state_class = finalize_init(init_steady_state_class, tile); %assigns T_first_cell and potentially start_depth_steady_state
            end
            
            CURRENT = tile.TOP_CLASS;
            CURRENT.STATVAR.top_depth_rel2groundSurface = 0; %set initial surface to zero
            
            CURRENT.STATVAR.T = CURRENT.STATVAR.layerThick .*0 + tile.PARA.T_first_cell;            
            CURRENT = convert_units(CURRENT, tile);
            CURRENT = finalize_init(CURRENT, tile);
            [CURRENT, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = ...
                init_T_steady_state_TOP_CLASS(CURRENT, tile.PARA.T_first_cell, tile.PARA.start_depth_steady_state, tile.FORCING.PARA.heatFlux_lb);

            CURRENT.PARA.target_grid = tile.GRID.STATVAR.GRID;
            CURRENT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
            while ~isequal(CURRENT.NEXT, tile.BOTTOM_CLASS.NEXT)
                CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
                
                CURRENT.NEXT.STATVAR.T = CURRENT.NEXT.STATVAR.layerThick .*0 + tile.PARA.T_first_cell;            
                CURRENT.NEXT = convert_units(CURRENT.NEXT, tile);
                CURRENT.NEXT = finalize_init(CURRENT.NEXT, tile);
                [CURRENT.NEXT, T_end, thermCond_end, layerThick_end, start_depth_steady_state] = ...
                    init_T_steady_state(CURRENT.NEXT, T_end, start_depth_steady_state, thermCond_end, layerThick_end, tile.FORCING.PARA.heatFlux_lb);
                
                CURRENT.NEXT.PARA.target_grid = tile.GRID.STATVAR.GRID;
                CURRENT.NEXT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
                
                CURRENT = CURRENT.NEXT;
            end

            %7. assign interaction classes 
            CURRENT = tile.TOP_CLASS;
            
            while ~isequal(CURRENT.NEXT, tile.BOTTOM)
                ia_class = get_IA_class(class(CURRENT), class(CURRENT.NEXT));
                CURRENT.IA_NEXT = ia_class;
                CURRENT.IA_NEXT.PREVIOUS = CURRENT;
                CURRENT.IA_NEXT.NEXT = CURRENT.NEXT;
                CURRENT.NEXT.IA_PREVIOUS = CURRENT.IA_NEXT;
                
                finalize_init(CURRENT.IA_NEXT, tile); %Added Sebastian, defined in IA_BASE as empty
                
                CURRENT = CURRENT.NEXT;
            end
            
            %8. assign SNOW class
            snow_class_name = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_name;
            snow_class_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_index;
            
            if ~isempty(snow_class_name) && sum(isnan(snow_class_name))==0
                snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
                snow_class = snow_class{snow_class_index,1};
                
                tile.STORE.SNOW = copy(snow_class);
                %no convert_units for SNOW classes needed at this point!
                tile.STORE.SNOW = finalize_init(tile.STORE.SNOW, tile); %make this dependent on TILE!
            end
            
            %9. assign sleeping classes
            sleeping_classes = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_name;
            sleeping_classes_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_index; 

            for i=1:size(sleeping_classes,1)
                sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
                sc = sc{sleeping_classes_index(i,1),1};
                tile.STORE.SLEEPING{i,1} = copy(sc);
                tile.STORE.SLEEPING{i,1} = convert_units(tile.STORE.SLEEPING{i,1}, tile);
                tile.STORE.SLEEPING{i,1} = finalize_init(tile.STORE.SLEEPING{i,1}, tile);
                tile.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
            end
            
            %10. assign time, etc.
            tile.t = tile.FORCING.PARA.start_time;
            
            %11. assign LATERAL classes 
            tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
            tile.LATERAL = finalize_init(tile.LATERAL, tile);
            
            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            tile.OUT = finalize_init(tile.OUT, tile);
            
        end
        
        function tile = build_tile_update_forcing_out(tile)
%             PARA2 = tile.PARA;
%             fields = fieldnames(tile.RUN_INFO.TILE);
%             for i=1:size(fields,1)
%                 if ~strcmp(fields{i,1}, 'OUT') || ~strcmp(fields{i,1}, 'FORCING') 
%                     tile.(fields{i,1}) = tile.RUN_INFO.TILE.(fields{i,1});
%                 end
%             end
%             fields = fieldnames(PARA2);
%             for i=1:size(fields,1)
%                 tile.PARA.(fields{i,1}) = PARA2.(fields{i,1});
%             end
%             
%             %1. forcing
%             tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});
%             tile.FORCING = finalize_init(tile.FORCING, tile);
% 
%             %10. assign time, etc.
%             tile.t = tile.FORCING.PARA.start_time;
% 
%             %12. assign OUT classes
%             tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
%             tile.OUT = finalize_init(tile.OUT, tile);
            
            %2. forcing -> special forcing class required
            tile.FORCING = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.forcing_class){tile.PARA.forcing_class_index,1});

            %12. assign OUT classes
            tile.OUT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.out_class){tile.PARA.out_class_index,1});
            
            %use old tile
            tile.GRID = tile.RUN_INFO.TILE.GRID;   
            tile.BOTTOM = tile.RUN_INFO.TILE.BOTTOM;
            tile.TOP = tile.RUN_INFO.TILE.TOP;
            tile.BOTTOM_CLASS = tile.RUN_INFO.TILE.BOTTOM_CLASS;
            tile.TOP_CLASS = tile.RUN_INFO.TILE.TOP_CLASS;
            tile.timestep = tile.RUN_INFO.TILE.timestep;
            tile.LATERAL = tile.RUN_INFO.TILE.LATERAL;     
            tile.STORE = tile.RUN_INFO.TILE.STORE;     
            
            %use old PARA, but overwrite all newly set values
            PARA_new = tile.PARA;
            tile.PARA = tile.RUN_INFO.TILE.PARA;
            fn = fieldnames(PARA_new);
            for i=1:size(fn,1)  %be careful, does not work if empty array (and not NaN) is willingly assigned to a parameter
                if ~isempty(PARA_new.(fn{i,1}))
                    tile.PARA.(fn{i,1}) = PARA_new.(fn{i,1});
                end
            end

            
            tile.FORCING = finalize_init(tile.FORCING, tile); 
            tile.OUT = finalize_init(tile.OUT, tile);           
            %10. assign time, etc.
            tile.TEMP.time_difference = tile.RUN_INFO.TILE.t - tile.FORCING.PARA.start_time; %Used to correct time variables in subsurface classes 
            tile.t = tile.FORCING.PARA.start_time;
            
            %reset IA time
            tile.LATERAL.IA_TIME = tile.t + tile.LATERAL.IA_TIME_INCREMENT;
            
            %reset time for BGC class (do mothing if no BGC class exists)
            %-> MAKE THIS A GENERAL RESET_TIME OR ADJUST_TIME FUNCTION THAT
            %IS DEFINED IN BASE AND OVERWRITTEN IN ALL FUNCTIONS THAT
            %ACTUALLY HAVE A TIME VARIABLE - CALCULATE TIME OFFSET BETWEEN
            %OLD AND NEW FORCING, I.E LAST TIMESTAMP OF OLD RUN AND FIRST TIMESTAMP OF NEW RUN 
            CURRENT = tile.TOP.NEXT;
            while ~isequal(CURRENT.NEXT, tile.BOTTOM)
                CURRENT = reset_timestamps(CURRENT, tile);
                CURRENT = CURRENT.NEXT;
            end
            
%             tile.LATERAL = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.lateral_class){tile.PARA.lateral_class_index,1});
%             tile.LATERAL = finalize_init(tile.LATERAL, tile);
            
            tile.RUN_INFO.TILE = tile;

            
            tile.PARA.run_name =  tile.RUN_INFO.PPROVIDER.PARA.run_name;
            tile.PARA.result_path =  tile.RUN_INFO.PPROVIDER.PARA.result_path;
        end
        
        
        function tile = build_tile_restart_OUT_last_timestep(tile)
            temp=load([tile.PARA.restart_file_path tile.PARA.restart_file_name]);
            variables = fieldnames(temp.out.STRATIGRAPHY);
            for i=1:size(variables,1)
                tile.(variables{i,1}) = temp.out.STRATIGRAPHY.(variables{i,1});
            end
%             tile.LATERAL.IA_CLASSES = {};
%             tile.LATERAL.PARA.num_realizations = 1;
%             tile.LATERAL.PARA.worker_number = 1;
%             tile.OUT.OUTPUT_TIME = tile.OUT.OUTPUT_TIME+100;
            %tile.LATERAL.IA_TIME = tile.FORCING.PARA.end_time;
        end
        
        
        
        %-------------param file generation-----
       function tile = param_file_info(varargin)

            if nargin==2
                tile = varargin{1};
                option = varargin{2};
            else
                tile = varargin{1};
                option = 0;
            end
            tile.PARA.STATVAR = [];
            tile.PARA.class_category = 'TILE';
            tile.PARA.default_value=[];
            tile.PARA.options=[];
            tile.PARA.comment=[];
            
            if strcmp(option, 'new_init')
                tile.PARA.builder = [];
                tile.PARA.default_value.builder = {'new_init'};
                parameters = {'latitude'; 'longitude'; 'altitude'; 'domain_depth'; 'area'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
                    'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
                    'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'};
                
                for i=1:size(parameters,1)
                    tile.PARA.(parameters{i,1})=[];
                end
                
                tile.PARA.comment.latitude = {'geographic coordinate, e.g. 70.956'};
                tile.PARA.comment.longitude = {'geographic coordinate, e.g. -158.123'};
                tile.PARA.comment.altitude = {'altitude [m]'};
                tile.PARA.default_value.domain_depth = {100};
                tile.PARA.comment.domain_depth = {'vertical depth of the model domain [m]'};
                tile.PARA.default_value.area = {1};
                tile.PARA.comment.area = {'area of the model domain [m2]'};
                tile.PARA.default_value.forcing_class = {'FORCING_seb'};
                tile.PARA.default_value.forcing_class_index = {1};
                tile.PARA.default_value.grid_class = {'GRID_user_defined'};
                tile.PARA.default_value.grid_class_index = {1};
                tile.PARA.default_value.out_class = {'OUT_all_lateral'};
                tile.PARA.default_value.out_class_index = {1};
                tile.PARA.default_value.strat_classes_class = {'STRAT_classes'};
                tile.PARA.default_value.strat_classes_class_index = {1};
                
                tile.PARA.comment.strat_statvar_class = {'list of STRATIGRAPHY_STATVAR classes that provide initial state of state variables'};
                tile.PARA.options.strat_statvar_class.name =  'H_LIST'; %
                tile.PARA.options.strat_statvar_class.entries_x = {'STRAT_layers' 'STRAT_linear'};
                tile.PARA.options.strat_statvar_class_index.name =  'H_LIST'; 
                tile.PARA.options.strat_statvar_class_index.entries_x = {1 1};
                
                tile.PARA.comment.lateral_class = {'lateral class, e.g. LATERAL1D or LATERAL3D'};
                tile.PARA.default_value.lateral_class = {'LATERAL_1D'};
                tile.PARA.default_value.lateral_class_index = {1};
                
                tile.PARA.comment.lateral_IA_classes = {'list of lateral interaction classes'};
                tile.PARA.options.lateral_IA_classes.name =  'H_LIST'; %
                tile.PARA.options.lateral_IA_classes_index.name =  'H_LIST';
                
            elseif strcmp(option, 'new_init_steady_state')
                tile.PARA.builder = [];
                tile.PARA.default_value.builder = {'new_init_steady_state'};
                parameters = {'latitude'; 'longitude'; 'altitude'; 'domain_depth'; 'area'; 'forcing_class'; 'forcing_class_index'; 'grid_class'; 'grid_class_index'; 'out_class'; ...
                    'out_class_index'; 'strat_classes_class'; 'strat_classes_class_index'; 'strat_statvar_class'; 'strat_statvar_class_index'; 'lateral_class'; ...
                    'lateral_class_index'; 'lateral_IA_classes'; 'lateral_IA_classes_index'; 'init_steady_state_class'; 'init_steady_state_class_index';...
                    'T_first_cell'; 'start_depth_steady_state'};
                for i=1:size(parameters,1)
                    tile.PARA.(parameters{i,1})=[];
                end
                
                tile.PARA.comment.latitude = {'geographic coordinate, e.g. 70.956'};
                tile.PARA.comment.longitude = {'geographic coordinate, e.g. -158.123'};
                tile.PARA.comment.altitude = {'altitude [m]'};
                tile.PARA.default_value.domain_depth = {100};
                tile.PARA.comment.domain_depth = {'vertical depth of the model domain [m]'};
                tile.PARA.default_value.area = {1};
                tile.PARA.comment.area = {'area of the model domain [m2]'};
                tile.PARA.default_value.forcing_class = {'FORCING_seb'};
                tile.PARA.default_value.forcing_class_index = {1};
                tile.PARA.default_value.grid_class = {'GRID_user_defined'};
                tile.PARA.default_value.grid_class_index = {1};
                tile.PARA.default_value.out_class = {'OUT_all_lateral'};
                tile.PARA.default_value.out_class_index = {1};
                tile.PARA.default_value.strat_classes_class = {'STRAT_classes'};
                tile.PARA.default_value.strat_classes_class_index = {1};
                
                tile.PARA.comment.strat_statvar_class = {'list of STRATIGRAPHY_STATVAR classes that provide initial state of state variables'};
                tile.PARA.options.strat_statvar_class.name =  'H_LIST'; %
                tile.PARA.options.strat_statvar_class.entries_x = {'STRAT_layers'};
                tile.PARA.options.strat_statvar_class_index.name =  'H_LIST'; 
                tile.PARA.options.strat_statvar_class_index.entries_x = {1};
                
                tile.PARA.comment.lateral_class = {'lateral class, e.g. LATERAL1D or LATERAL3D'};
                tile.PARA.default_value.lateral_class = {'LATERAL_1D'};
                tile.PARA.default_value.lateral_class_index = {1};
                
                tile.PARA.comment.lateral_IA_classes = {'list of lateral interaction classes'};
                tile.PARA.options.lateral_IA_classes.name =  'H_LIST'; %
                tile.PARA.options.lateral_IA_classes_index.name =  'H_LIST';
                
                tile.PARA.comment.init_steady_state_class = {'init_steady_state class to compute temperature of first grid cell, leave empty when using T_first_grid_cell'};
                tile.PARA.comment.T_first_cell = {'temperature of first grid cell used to compute the temperature gradient, leave empty when using init_steady_state class'};
                tile.PARA.comment.start_depth_steady_state = {'depth [m] where temperature gradient starts, constant above, leave empty when using init_steady_state class'};
                
            elseif strcmp(option, 'update_forcing_out')
                tile.PARA.builder = [];
                tile.PARA.default_value.builder = {'update_forcing_out'};
                parameters = { 'forcing_class'; 'forcing_class_index';  'out_class'; 'out_class_index'};
                
                for i=1:size(parameters,1)
                    tile.PARA.(parameters{i,1})=[];
                end
                
                tile.PARA.default_value.forcing_class = {'FORCING_seb'};
                tile.PARA.default_value.forcing_class_index = {1};
                tile.PARA.default_value.out_class = {'OUT_all_lateral'};
                tile.PARA.default_value.out_class_index = {1};
                
            elseif strcmp(option, 'restart_OUT_last_timestep')
                tile.PARA.builder = [];
                tile.PARA.default_value.builder = {'restart_OUT_last_timestep'};
                parameters = {'restart_file_path'; 'restart_file_name'};
                for i=1:size(parameters,1)
                    tile.PARA.(parameters{i,1})=[];
                end
                
                tile.PARA.comment.restart_file_path = {'path and filename of restart file'};
            else
                tile = provide_PARA(tile);
            end

        end
        

    end
end



