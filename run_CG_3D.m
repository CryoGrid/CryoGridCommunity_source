%========================================================================
% CryoGrid main run file, one-dimensional (single-tile) configuration
% S. Westermann, T. Ingemann-Nielsen, J. Scheer, October 2020
%========================================================================

clear all

%=========================================================================
%MODIFY TO RUN
%=========================================================================
%set parameter files for initialization
init_format = 'EXCEL'; %EXCEL or YAML
run_number = 'Herschell'; %paramter file name and result directory 
const_file = 'CONSTANTS_excel'; %file with constants
%result_path = '../CryoGrid_Git_results/';
result_path = '../results/';

%must be part of the overarching concept to initialize multi-tile runs 
lateral.PARA.run_number = [1;2;3];
number_of_tiles = 3;

%=========================================================================
%DO NOT MODIFY BELOW
%=========================================================================
%=========================================================================
%INITIALIZATION
%=========================================================================
%set various paths and start parallel pool
modules_path = 'modules';
addpath(genpath(modules_path));
config_path = fullfile(result_path, run_number);
forcing_path = fullfile ('./forcing/');
parpool(number_of_tiles)



spmd
    
    %set run number and parameter file name for own worker - this should become part of the overarching concept to initialize multi-tile runs 
    run_number = [run_number '_' num2str(lateral.PARA.run_number(labindex,1))]; 
    parameter_file = [run_number '.xlsx'];
    
    %call provider classes to extract run parameters and constants from the
    %parameter/constant files
    pprovider = PARAMETER_PROVIDER(config_path, init_format, run_number);
    cprovider = CONSTANT_PROVIDER(config_path, init_format, const_file);
    fprovider = FORCING_PROVIDER(pprovider, forcing_path);

    % Build the model tile (forcing, grid, out and stratigraphy classes)
    %tile_built = TILE_BUILDER(pprovider, cprovider, fprovider);
    tile = TILE(pprovider, cprovider, fprovider, run_number, result_path);
    tile.LATERAL = LATERAL_3D(tile);
    
%     tile = TILE();
%     tile.CONST.day_sec = 24.*3600;

    
    %initialize LATERAL classes as defined in the parameter file

    %lateral = LATERAL_3D(tile_built);

    %this should be cleaned up and become part of the initialization
    %procedure for multi-tile runs
    tile.LATERAL = assign_number_of_realizations(tile.LATERAL, number_of_tiles);
    tile.LATERAL = get3d_PARA(tile.LATERAL);
    tile.LATERAL = get_index(tile.LATERAL);

    %assign forcing and out classes
%     tile.FORCING = tile_built.forcing;
%     tile.OUT = tile_built.out;
% 
%     tile.RUN_NUMBER = run_number;
%     tile.RESULT_PATH = result_path;
    
    %assign bottom and top classes
    TOP_CLASS = tile.TOP_CLASS;
    TOP = tile.TOP;
    %tile.TOP = TOP;

    BOTTOM_CLASS = tile.BOTTOM_CLASS;
    BOTTOM = tile.BOTTOM;
    %tile.BOTTOM = BOTTOM;

    TOP.LATERAL = tile.LATERAL;
%     TOP.FORCING = forcing;
%     TOP.RUN_NUMBER = run_number;
%     TOP.RESULT_PATH = result_path;

    %global variable assigned here (possibly remove this in the future)
%     day_sec = 24.*3600;
    
    %initialize running time variable t [days]
    %tile.t = tile.FORCING.PARA.start_time;
    
    %lateral = initialize_lateral_3D(lateral, TOP, BOTTOM, t);

  %  clear pprovider cprovider fprovider forcing_path config_path modules_path run_number result_path const_file init_format
    %lateral.IA_TIME = t;

    %lkjlkjlkj
    
    while tile.t < tile.FORCING.PARA.end_time
        
        tile = interpolate_forcing_tile(tile);
        %---------boundary conditions
        
        %proprietary function for each class, i.e. the "real upper boundary"
        %only evaluated for the first cell/block
        
        TOP.NEXT = get_boundary_condition_u(TOP.NEXT, tile);
        CURRENT = TOP.NEXT;
        
        %function independent of classes, each class must comply with this function!!!
        %evaluated for every interface between two cells/blocks
        while ~isequal(CURRENT.NEXT, BOTTOM)
            get_boundary_condition_m(CURRENT.IA_NEXT, tile);
            CURRENT = CURRENT.NEXT;
        end
        
        %proprietary function for each class, i.e. the "real lower boundary"
        %only evaluated for the last cell/block
        CURRENT = get_boundary_condition_l(CURRENT, tile);  %At this point, CURRENT is equal to BOTTOM_CLASS
        %--------------------------
        
        %calculate spatial derivatives for every cell in the stratigraphy
        CURRENT = TOP.NEXT;
        while ~isequal(CURRENT, BOTTOM)
            CURRENT = get_derivatives_prognostic(CURRENT, tile);
            CURRENT = CURRENT.NEXT;
        end
        
        %calculate minimum timestep required for all cells in days
        CURRENT = TOP.NEXT;
        tile.timestep = tile.CONST.day_sec;
        while ~isequal(CURRENT, BOTTOM)
            tile.timestep = min(tile.timestep, get_timestep(CURRENT, tile));
            CURRENT = CURRENT.NEXT;
        end
        
        tile.next_break_time = min(tile.LATERAL.IA_TIME, tile.OUT.OUTPUT_TIME);
        tile.timestep = min(tile.timestep, (tile.next_break_time - tile.t).*tile.CONST.day_sec);
        %make sure to hit the output times!
        
        %calculate prognostic variables
        CURRENT = TOP.NEXT;
        while ~isequal(CURRENT, BOTTOM)
            CURRENT = advance_prognostic(CURRENT, tile);
            CURRENT = CURRENT.NEXT;
        end
                
        %calculate diagnostic variables
        %some effects only happen in the first cell
        TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, tile);
        
        CURRENT = BOTTOM.PREVIOUS;
        while ~isequal(CURRENT, TOP)
            CURRENT = compute_diagnostic(CURRENT, tile);
            CURRENT = CURRENT.PREVIOUS;
        end
        
        %check for triggers that reorganize the stratigraphy
        CURRENT = TOP.NEXT;
        while ~isequal(CURRENT, BOTTOM)
            CURRENT = check_trigger(CURRENT, tile);
            CURRENT = CURRENT.NEXT;
        end

        TOP_CLASS = TOP.NEXT; %TOP_CLASS and BOTTOM_CLASS for convenient access
        BOTTOM_CLASS = BOTTOM.PREVIOUS;
        %TOP.TIME = tile.t;
        
        tile = interact_lateral(tile);
        %lateral = interact(lateral, forcing, t);
        %lateral = lateral_IA(lateral, forcing, t);

        %calculate new time
        tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
        
        %store the output according to the defined OUT clas
        tile = store_OUT_tile(tile);        
        %out = store_OUT(out, t, TOP, BOTTOM, timestep);        
        
    end
    
end

delete(gcp('nocreate'));


