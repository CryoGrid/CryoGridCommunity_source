%========================================================================
% CryoGrid main run file, one-dimensional (single-tile) configuration
% S. Westermann, T. Ingemann-Nielsen, J. Scheer, October 2020
%========================================================================

clear all

%=========================================================================
%MODIFY TO RUN
%=========================================================================
% %set parameter files for initialization
% init_format = 'EXCEL'; %EXCEL or YAML
% run_number = 'example4'; %paramter file name and result directory 
% const_file = 'CONSTANTS_excel'; %file with constants
% %result_path = '../CryoGrid_Git_results/';
% result_path = '../results/';

%set parameter files for initialization
init_format = 'EXCEL'; %EXCEL or YAML
run_number = 'example2'; %paramter file name and result directory 
%run_number = 'example1'; %paramter file name and result directory 
const_file = 'CONSTANTS_excel'; %file with constants
%result_path = '../CryoGrid_Git_results/';
result_path = './results/';

%=========================================================================
%DO NOT MODIFY BELOW
%=========================================================================
%=========================================================================
%INITIALIZATION
%=========================================================================
%set various paths
modules_path = 'modules';
addpath(genpath(modules_path));
config_path = fullfile(result_path, run_number);
forcing_path = fullfile ('./forcing/');

%call provider classes to extract run parameters and constants from the
%parameter/constant files
pprovider = PARAMETER_PROVIDER(config_path, init_format, run_number);
cprovider = CONSTANT_PROVIDER(config_path, init_format, const_file);
fprovider = FORCING_PROVIDER(pprovider, forcing_path);

% Build the model tile (forcing, grid, out and stratigraphy classes)
tile = TILE(pprovider, cprovider, fprovider, run_number, result_path);

%assign bottom and top classes
TOP_CLASS = tile.TOP_CLASS;
BOTTOM_CLASS = tile.BOTTOM_CLASS;
TOP = tile.TOP;
BOTTOM = tile.BOTTOM;
TOP.LATERAL = tile.LATERAL;

clear pprovider cprovider fprovider forcing_path config_path modules_path run_number result_path const_file init_format

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
    %lateral = interact(lateral, forcing, t);
    tile = interact_lateral(tile);
    
    %set TOP_CLASS and BOTTOM_CLASS for convenient access
    TOP_CLASS = TOP.NEXT; 
    BOTTOM_CLASS = BOTTOM.PREVIOUS;

    
    %update time variable t
    tile.t = tile.t + tile.timestep./tile.CONST.day_sec;
    
    %model 
    %out = store_OUT(out, t, TOP, BOTTOM, timestep);
    tile = store_OUT_tile(tile);
end

