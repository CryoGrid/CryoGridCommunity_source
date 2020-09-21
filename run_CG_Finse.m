% CryoGrid main file to be executed
% S.Westermann, Jan 2019


clear all
modules_path = 'modules';
addpath(genpath(modules_path));

run_number = 'Finse_4432';      
result_path = './results/';
config_path = fullfile(result_path, run_number);
forcing_path = fullfile ('./forcing/');

parameter_file = [run_number '.xlsx'];
const_file = 'CONSTANTS_excel.xlsx';              

% =====================================================================
% Use modular interface to build model run
% =====================================================================

% Depending on parameter_file_type, instantiates
% PARAMETER_PROVIDER, CONSTANT_PROVIDER and FORCING_PROVIDER classes

pprovider = PARAMETER_PROVIDER_EXCEL(config_path, parameter_file);
cprovider = CONSTANT_PROVIDER_EXCEL(config_path, const_file);
fprovider = FORCING_PROVIDER(pprovider, forcing_path);

% Build the actual model tile (forcing, grid, out and stratigraphy classes)
tile = TILE_BUILDER(pprovider, cprovider, fprovider);
lateral = LATERAL_1D(tile);

% ------ preparations ----------------------

forcing = tile.forcing;
out = tile.out;

TOP_CLASS = tile.TOP_CLASS;
BOTTOM_CLASS = tile.BOTTOM_CLASS;
TOP = tile.TOP;
BOTTOM = tile.BOTTOM;

day_sec = 24.*3600;
t = forcing.PARA.start_time;
%t is in days, timestep should also be in days

% ------ time integration ------------------

%t is in days, timestep should also be in days

while t < forcing.PARA.end_time
    
    forcing = interpolate_forcing(t, forcing);
    %---------boundary conditions
    
    %proprietary function for each class, i.e. the "real upper boundary"
    %only evaluated for the first cell/block
    
    TOP.NEXT = get_boundary_condition_u(TOP.NEXT, forcing);
    CURRENT = TOP.NEXT;
   
    %function independent of classes, each class must comply with this function!!!
    %evaluated for every interface between two cells/blocks
    while ~isequal(CURRENT.NEXT, BOTTOM)
        get_boundary_condition_m(CURRENT.IA_NEXT);
        CURRENT = CURRENT.NEXT;
    end
    %proprietary function for each class, i.e. the "real lower boundary"
    %only evaluated for the last cell/block
    CURRENT = get_boundary_condition_l(CURRENT,  forcing);  %At this point, CURRENT is equal to BOTTOM_CLASS
    %--------------------------
    
    %calculate spatial derivatives for every cell in the stratigraphy
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = get_derivatives_prognostic(CURRENT);
        CURRENT = CURRENT.NEXT;
    end
    
    %calculate minimum timestep required for all cells in days
    CURRENT = TOP.NEXT;
    timestep=3600;
    while ~isequal(CURRENT, BOTTOM)
        timestep = min(timestep, get_timestep(CURRENT));
        CURRENT = CURRENT.NEXT;
    end
    next_break_time = min(lateral.IA_TIME, out.OUTPUT_TIME);
    timestep = min(timestep, (next_break_time - t).*day_sec);
    %make sure to hit the output times!
    
    %calculate prognostic variables
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = advance_prognostic(CURRENT, timestep);
        
    CURRENT = CURRENT.NEXT;
    end
    
    
    %calculate diagnostic variables
    %some effects only happen in the first cell
    TOP.NEXT = compute_diagnostic_first_cell(TOP.NEXT, forcing);
    if isnan(TOP.NEXT.STATVAR.Lstar)
        keyboard
    end
    
    CURRENT = BOTTOM.PREVIOUS;
    while ~isequal(CURRENT, TOP)
        CURRENT = compute_diagnostic(CURRENT, forcing);
        CURRENT = CURRENT.PREVIOUS;
    end
    
    %check for triggers that reorganize the stratigraphy
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = check_trigger(CURRENT, forcing);
        CURRENT = CURRENT.NEXT;
    end
    
    lateral = interact(lateral, forcing, t);
    
    
    TOP_CLASS = TOP.NEXT; %TOP_CLASS and BOTOOM_CLASS for convenient access
    BOTTOM_CLASS = BOTTOM.PREVIOUS;
    
       
    %calculate new time
    t = t + timestep./day_sec;
    
    %store the output according to the defined OUT clas
    out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path);
    
end

