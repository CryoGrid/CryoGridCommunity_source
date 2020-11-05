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
run_number = 'example1'; %paramter file name and result directory 
const_file = 'CONSTANTS_excel'; %file with constants
result_path = '../CryoGrid_Git_results/';
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
tile = TILE_BUILDER(pprovider, cprovider, fprovider);

%initialize LATERAL classes as defined in the parameter file
lateral = LATERAL_1D(tile);

%assign focing and out classes
forcing = tile.forcing;
out = tile.out;

%assign bottom and top classes
TOP_CLASS = tile.TOP_CLASS;
BOTTOM_CLASS = tile.BOTTOM_CLASS;
TOP = tile.TOP;
BOTTOM = tile.BOTTOM;    
TOP.LATERAL = lateral;

%global variable assigned here (possibly remove this in the future) 
day_sec = 24.*3600;

%initialize running time variable t [days]
t = forcing.PARA.start_time;

%=========================================================================
%TIME INTEGRATION
%=========================================================================
while t < forcing.PARA.end_time
    
    %interpolate focing data to time t
    forcing = interpolate_forcing(t, forcing);

    %upper boundar condition (uppermost class only)
    TOP.NEXT = get_boundary_condition_u(TOP.NEXT, forcing);
    
    %set fluxes between classes in the stratigrapht
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT.NEXT, BOTTOM)
        get_boundary_condition_m(CURRENT.IA_NEXT); %call interaction class function
        CURRENT = CURRENT.NEXT;
    end

    %lower boundary condition (lowermost class)
    CURRENT = get_boundary_condition_l(CURRENT,  forcing);  %At this point, CURRENT is equal to BOTTOM_CLASS

    %calculate spatial derivatives
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = get_derivatives_prognostic(CURRENT);
        CURRENT = CURRENT.NEXT;
    end
    
    %calculate timestep [second]
    CURRENT = TOP.NEXT;
    timestep = day_sec;
    while ~isequal(CURRENT, BOTTOM)
        timestep = min(timestep, get_timestep(CURRENT));
        CURRENT = CURRENT.NEXT;
    end
    next_break_time = min(lateral.IA_TIME, out.OUTPUT_TIME);
    timestep = min(timestep, (next_break_time - t).*day_sec);
    
    %prognostic step - integrate prognostic variables in time
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = advance_prognostic(CURRENT, timestep);
        CURRENT = CURRENT.NEXT;
    end
           
    %diagnostic step - compute diagnostic variables
    TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, forcing); %calculate Lstar, only uppermost class
    CURRENT = BOTTOM.PREVIOUS;
    while ~isequal(CURRENT, TOP)
        CURRENT = compute_diagnostic(CURRENT, forcing);
        CURRENT = CURRENT.PREVIOUS;
    end
    
    %triggers
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = check_trigger(CURRENT, forcing);
        CURRENT = CURRENT.NEXT;
    end
    
    %lateral interactions
    lateral = interact(lateral, forcing, t);
    
    %set TOP_CLASS and BOTTOM_CLASS for convenient access
    TOP_CLASS = TOP.NEXT; 
    BOTTOM_CLASS = BOTTOM.PREVIOUS;
    
    %update time variable t
    t = t + timestep./day_sec;
    
    %model output
    out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path);
end

