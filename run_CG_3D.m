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
result_path = './results/';
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
    run_number = [run_number '_' num2str(lateral.PARA.run_number(labidex,1))]; 
    parameter_file = [run_number '.xlsx'];

    %call provider classes to extract run parameters and constants from the
    %parameter/constant files
    pprovider = PARAMETER_PROVIDER_EXCEL(config_path, parameter_file);
    cprovider = CONSTANT_PROVIDER_EXCEL(config_path, const_file);
    fprovider = FORCING_PROVIDER(forcing_path, forcing_file);

    % Build the model tile (forcing, grid, out and stratigraphy classes)
    tile = TILE_BUILDER(pprovider, cprovider, fprovider);
    
    %initialize LATERAL classes as defined in the parameter file
    lateral = LATERAL_3D(tile);
    %this should be cleaned up and become part of the initialization
    %procedure for multi-tile runs
    lateral = assign_number_of_realizations(lateral, number_of_tiles);
    lateral = get3d_PARA(lateral);
    lateral = get_index(lateral);
    
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
    
    %lateral = initialize_lateral_3D(lateral, TOP, BOTTOM, t);

    
    %lateral.IA_TIME = t;

    %lkjlkjlkj
    
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
        
%         if t> datenum(1996,5,9,18,0,0) && lateral.STATVAR.index == 1
%             disp('Hallo1')
%         end
        
        
        CURRENT = BOTTOM.PREVIOUS;
        while ~isequal(CURRENT, TOP)
            CURRENT = compute_diagnostic(CURRENT, forcing);
            CURRENT = CURRENT.PREVIOUS;
        end
%         
%         if t> datenum(1996,5,9,18,0,0) && lateral.STATVAR.index == 1
%             disp('Hallo2')
%         end
        
        %check for triggers that reorganize the stratigraphy
        CURRENT = TOP.NEXT;
        while ~isequal(CURRENT, BOTTOM)
            CURRENT = check_trigger(CURRENT, forcing);
            CURRENT = CURRENT.NEXT;
        end
% 
%         if t> datenum(1996,5,9,18,0,0) && lateral.STATVAR.index == 1
%             disp('Hallo3')
%         end
        
        TOP_CLASS = TOP.NEXT; %TOP_CLASS and BOTTOM_CLASS for convenient access
        BOTTOM_CLASS = BOTTOM.PREVIOUS;
        
        
        %calculate new time
        t = t + timestep./day_sec;
        
        lateral = lateral_IA(lateral, forcing, t);
        
%         if t> datenum(1996,5,9,18,0,0) && lateral.STATVAR.index == 1
%             disp('Hallo4')
%         end
        
        %store the output according to the defined OUT clas
        out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number, timestep, result_path);
        
%         if t> datenum(1996,5,9,18,0,0) && lateral.STATVAR.index == 1
%             disp('Hallo5')
%         end
        
    end
    
end

delete(gcp('nocreate'));


