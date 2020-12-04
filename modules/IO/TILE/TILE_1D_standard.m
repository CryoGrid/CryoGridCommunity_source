% base class build a model tile

classdef TILE_1D_standard < matlab.mixin.Copyable
    
    properties
        
%         RUN_NUMBER
%         RESULT_PATH
        
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
        TOP          % Indicates current top position
        BOTTOM       % Indicates current bottom position
        TOP_CLASS    % Uppermost ground/snow class in the soil column
        BOTTOM_CLASS % Lowermost ground class in the soil column

    end
    
    
    methods
        
       

        function tile = provide_PARA(tile)

            tile.PARA.coordinates = [];
            tile.PARA.crs = [];
            tile.PARA.height_system = [];
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
            
        end
        
        function tile = provide_CONST(tile)

            tile.CONST.day_sec = [];

        end
        
        function tile = provide_STATVAR(tile)

        end
        

        function tile = initialize_excel(tile)
            
        end
       
        
        %assemble the stratigraphy
        function tile = finalize_init(tile)
            
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
            class_depths = [class_depths; tile.GRID.GRID(end,1)];
            
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
            CURRENT = finalize_init(CURRENT, tile);

            CURRENT.PARA.target_grid = tile.GRID.GRID;
            CURRENT.PARA.target_layerThick = tile.GRID.STATVAR.layerThick;
            while ~isequal(CURRENT.NEXT, tile.BOTTOM_CLASS.NEXT)
                CURRENT.NEXT.STATVAR.top_depth_rel2groundSurface = CURRENT.STATVAR.top_depth_rel2groundSurface + sum(CURRENT.STATVAR.layerThick,1);
                
                CURRENT.NEXT = finalize_init(CURRENT.NEXT, tile);

                CURRENT.NEXT.PARA.target_grid = tile.GRID.GRID;
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
                
                
                CURRENT = CURRENT.NEXT;
            end
            
            %8. assign SNOW class
            snow_class_name = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_name;
            snow_class_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.snow_class_index;
            
            snow_class =  tile.RUN_INFO.PPROVIDER.CLASSES.(snow_class_name);
            snow_class = snow_class{snow_class_index,1}; 
           
            tile.TOP.STORE.SNOW = copy(snow_class);
            tile.TOP.STORE.SNOW = finalize_init(tile.TOP.STORE.SNOW, tile); %make this dependent on TILE!
            
            %9. assign sleeping classes
            sleeping_classes = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_name;
            sleeping_classes_index = tile.RUN_INFO.PPROVIDER.CLASSES.(tile.PARA.strat_classes_class){tile.PARA.strat_classes_class_index,1}.PARA.sleeping_classes_index; 

            for i=1:size(sleeping_classes,1)
                sc = tile.RUN_INFO.PPROVIDER.CLASSES.(sleeping_classes{i,1});
                sc = sc{sleeping_classes_index(i,1),1};
                tile.TOP.STORE.SLEEPING{i,1} = copy(sc);
                tile.TOP.STORE.SLEEPING{i,1} = finalize_init(tile.TOP.STORE.SLEEPING{i,1}, tile);
                tile.TOP.STORE.SLEEPING{i,2} = sleeping_classes_index(i,1);
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
        

        function tile = interpolate_forcing_tile(tile)
             tile.FORCING = interpolate_forcing(tile.t, tile.FORCING);
        end

        function tile = interact_lateral(tile)
            tile.LATERAL = interact(tile.LATERAL, tile);
        end
        
        function tile = store_OUT_tile(tile)
            tile.OUT = store_OUT(tile.OUT, tile);
        end        
        
        
        
        function tile = run(tile)
            
            TOP_CLASS = tile.TOP_CLASS;
            BOTTOM_CLASS = tile.BOTTOM_CLASS;
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
                
                tile = interact_lateral(tile);
                
                %set TOP_CLASS and BOTTOM_CLASS for convenient access
                TOP_CLASS = TOP.NEXT;
                BOTTOM_CLASS = BOTTOM.PREVIOUS;
                
                %update time variable t
                tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
                
                %model
                tile = store_OUT_tile(tile);
            end
            
        end
        

    end
end



