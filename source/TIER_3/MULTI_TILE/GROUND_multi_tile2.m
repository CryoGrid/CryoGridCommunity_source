%========================================================================
% CryoGrid GROUND class 

% S. Westermann, October 2020
%========================================================================

classdef GROUND_multi_tile2 < SEB  %put all required ground classes here

%     properties
%         SUB_TILES_TOP
%         SUB_TILES_BOTTOM
%     end
    
    methods
        
        %----mandatory functions---------------
        %----initialization--------------------
        
        
        function ground = provide_PARA(ground)
            ground.PARA.areal_fractions = []; %total of 1
            ground.PARA.surface_elevation = []; %offset from tile surface elevation
            ground.PARA.contact_lengths = []; %matrix
            ground.PARA.distances = []; %matrix
            ground.PARA.ground_classes = []; %vector
            ground.PARA.ground_classes_index = []; %vector
            ground.PARA.grid_classes = []; %vector (one grid class per ground class)
            ground.PARA.grid_classes_index = []; %vector (one grid class per ground class)
            ground.PARA.strat_statvar_classes = []; %vector (as many as needed)
            ground.PARA.strat_statvar_classes_index = []; %vector (as many as needed)
            ground.PARA.strat_statvar_classes_matrix = []; %matrix (statvar class x ground class)
            
            ground.PARA.unit_conversion_class = 'UNIT_CONVERSION_standard'; %can be overwritten if needed
        end
        
        function ground = provide_CONST(ground)
            ground.CONST.cp = [];
            ground.CONST.kappa = [];
            ground.CONST.g = [];
            ground.CONST.L_s = [];
        end
        
        function ground = provide_STATVAR(ground)
            ground.STATVAR.Lstar = -100;
            ground.STATVAR.Qh = 0;
            ground.STATVAR.Qe = 0;
        end
        
        function ground = convert_units(ground, tile)
            ground.PARA.airT_height = tile.FORCING.PARA.airT_height;
            %do nothing, handled in the subclasses
        end
        
        function ground = finalize_init(ground, tile)
            %2. grid
            ground.STATVAR.layerThick = 0;
            for ii = 1:size(ground.PARA.ground_classes,1)

                GRID = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ground.PARA.grid_classes{ii,1}){ground.PARA.grid_classes_index(ii,1),1});
                GRID = finalize_init_GROUND_multi_tile(GRID, tile);
                
                %3. map statvar to grid, using STRATIGRAPHY_STATVAR classes
                for i=1:size(ground.PARA.strat_statvar_classes_matrix,1)
                    if ground.PARA.strat_statvar_classes_matrix(i,ii)==1
                        strat_statvar_class = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ground.PARA.strat_statvar_classes{i,1}){ground.PARA.strat_statvar_classes_index(i,1),1});
                        strat_statvar_class = finalize_init_GROUND_multi_tile(strat_statvar_class, GRID);  %NEW FUNCTION NEEDED
                    end
                end
                
                %4. build stratigraphy
                ground.STATVAR.SUB_TILES_TOP{ii,1} = Top();
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT = copy(tile.RUN_INFO.PPROVIDER.CLASSES.(ground.PARA.ground_classes{ii,1}){ground.PARA.ground_classes_index(ii,1)});
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.PREVIOUS = ground.STATVAR.SUB_TILES_TOP{ii,1};
                ground.STATVAR.SUB_TILES_BOTTOM{ii,1} = Bottom();
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.NEXT = ground.STATVAR.SUB_TILES_BOTTOM{ii,1};
                ground.STATVAR.SUB_TILES_BOTTOM{ii,1}.PREVIOUS = ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT;
                
                %assigned to the main class before
                upperPos = ground.STATVAR.upperPos + ground.PARA.surface_elevation(ii,1);
                lowerPos  = ground.STATVAR.lowerPos;

                variables = fieldnames(ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.STATVAR);
                range = (GRID.STATVAR.MIDPOINTS <= upperPos-lowerPos);  % class_depths is upper and lower Pos
                for j=1:size(variables,1)
                    if isfield(GRID.STATVAR, variables{j,1})
                        ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.STATVAR.(variables{j,1}) = GRID.STATVAR.(variables{j,1})(range);
                    end
                end

                %6. set top depths relative to surface and finalize initialization for
                %subsurface classes
                total_area = tile.PARA.area;
                tile.PARA.area = tile.PARA.area .* ground.PARA.areal_fractions(ii,1); %adjust tarea
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT = convert_units(ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT, tile);
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT = finalize_init(ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT, tile);
                ground.STATVAR.layerThick = ground.STATVAR.layerThick + tile.PARA.area .* sum(ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.STATVAR.layerThick);
                
                tile.PARA.area = total_area;
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.PARA.target_grid = GRID.STATVAR.GRID;
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT.PARA.target_layerThick = GRID.STATVAR.layerThick;

                
            end
            ground.STATVAR.layerThick = ground.STATVAR.layerThick ./ tile.PARA.area; %average thickness of all subtile classes, needed in TILE
            %Depending on class, assign LUT based on single class 1  
            
        end
        
        
        %---time integration------
        
        function ground = get_boundary_condition_u(ground, tile)
            %ground = get_boundary_condition_u1(ground, tile);
            height_of_class = [];
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT;
                height = sum(CURRENT.STATVAR.layerThick,1);
                while ~strcmp(
            end
            
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT = get_boundary_condition_u(ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT, tile);
            end
            %ground = get_boundary_condition_u2(ground, tile);
        end
        
%         function [ground, S_up] = penetrate_SW(ground, S_down)  %mandatory function when used with class that features SW penetration
%             [ground, S_up] = penetrate_SW@GROUND_freezeC_bucketW_Xice_seb_snow(ground, S_down);
%         end
%NEED TO THINK THIS THROUGH! PROBABLY NOT NEEDED FOR MOST CASES
        
        
        function ground = get_boundary_condition_l(ground, tile) %function will generally not be calles
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
%                 CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT;
%                 while ~strcmp(class(CURRENT.NEXT), 'Bottom')
%                     CURRENT= CURRENT.NEXT;
%                 end
                CURRENT = get_boundary_condition_l(ground.STATVAR.SUB_TILES_BOTTOM{ii,1}.PREVIOUS, tile);
            end
        end
        
        
        function ground = get_derivatives_prognostic(ground, tile)
            %ground = get_derivatives_prognostic1(ground, tile);
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT;
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    get_boundary_condition_m(CURRENT.IA_NEXT, tile); %call interaction class function
                    CURRENT = CURRENT.NEXT;
                end
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1};
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    CURRENT= CURRENT.NEXT;
                    CURRENT = get_derivatives_prognostic(CURRENT, tile);
                end
            end
            %ground = get_derivatives_prognostic2(ground, tile);
        end
        
        function timestep = get_timestep(ground, tile)
            timestep = 1e8;
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1};
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    CURRENT= CURRENT.NEXT;
                    timestep = min(timestep, get_timestep(CURRENT, tile));
                end
            end
        end
        
        function ground = advance_prognostic(ground, tile)
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1};
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    CURRENT= CURRENT.NEXT;
                    CURRENT = advance_prognostic(CURRENT, tile);
                end
            end
        end
    
        
        function ground = compute_diagnostic_first_cell(ground, tile)
            %mix Qe and Qh here, and then compute bulk Obukhov length,
            %assign this to all sub tiles
            ground.STATVAR.Qe = 0;
            ground.STATVAR.Qh = 0;
            z0 = 0;
            area_acc = 0;
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1}.NEXT;
                ground.STATVAR.Qe = ground.STATVAR.Qe + CURRENT.STATVAR.Qh .* CURRENT.STATVAR.area(1,1);
                ground.STATVAR.Qh = ground.STATVAR.Qh + CURRENT.STATVAR.Qh .* CURRENT.STATVAR.area(1,1);
                z0 = z0 + CURRENT.PARA.z0 .* CURRENT.STATVAR.area(1,1);
                area_acc = area_acc + CURRENT.STATVAR.area(1,1);
            end
            ground.STATVAR.Qe = ground.STATVAR.Qe ./ area_acc;
            ground.STATVAR.Qh = ground.STATVAR.Qh ./ area_acc;
            ground.PARA.z0 = z0 ./ area_acc;
            ground = L_star(ground, tile.FORCING);
        end
        
        function ground = compute_diagnostic(ground, tile)
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1};
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    CURRENT= CURRENT.NEXT;
                    CURRENT = compute_diagnostic(CURRENT, tile);
                end
            end
        end
        
        function ground = check_trigger(ground, tile)
            for ii = 1:size(ground.STATVAR.SUB_TILES_TOP,1)
                CURRENT = ground.STATVAR.SUB_TILES_TOP{ii,1};
                while ~strcmp(class(CURRENT.NEXT), 'Bottom')
                    CURRENT= CURRENT.NEXT;
                    CURRENT = check_trigger(CURRENT, tile);
                end
            end
        end
        
    end
end