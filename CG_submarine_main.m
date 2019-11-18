
function ok = CG_submarine_main(run_number)

% CryoGrid main file to be executed
% based on S. Westermann, Jan 2019
%
% Adaptation for submarine permafrost
% F. Miesner Jan 2019

%clear all
addpath(genpath('modules'))


%results are stored in results/"run_number"
%this folder must contain an Excel spreadsheet called "run_number".xlsx deifining the properties of all classes requried for the run and the
%spreadsheet "CONSTANTS_excel.xlsx" containing all constants

if nargin < 1
    run_number = 'submarine_PF_test_salt'; %submarine_PF_test1
end
result_path = 'results/';
parameter_file_type = 'xlsx';
const_file = 'CONSTANTS_excel.xlsx';

%read information for all classes from Excel file and store it in cell
%arrays
parameter_info = read_excel2cell([result_path run_number '/' run_number '.' parameter_file_type]);  %read excel file to cell array

forcing = read_forcing_from_file(parameter_info);
grid = read_grid_from_file(parameter_info);
grid = reduce_grid(grid, forcing);
out = read_out_from_file(parameter_info);
out = complete_init_out(out, forcing);

stratigraphy_list = read_stratigraphies_from_file(parameter_info);
for i=1:size(stratigraphy_list,1)
    stratigraphy_list{i,1} = interpolate_to_grid(stratigraphy_list{i,1}, grid);
end

class_list = read_classes_from_file(parameter_info);
const_info = read_excel2cell([result_path run_number '/' const_file]);  %read CONST from excel file to cell array
for i=1:size(class_list,1)
    class_list{i,1}.CONST = initialize_from_file(class_list{i,1}, class_list{i,1}.CONST, const_info);
    class_list{i,1} = assign_global_variables(class_list{i,1}, forcing);
end
    
%assemble the model stratigraphy and define interactions between classes
[TOP_CLASS, BOTTOM_CLASS, TOP, BOTTOM] = assemble_stratigraphy(class_list, stratigraphy_list, grid, forcing);
TOP_CLASS = assemble_interactions(TOP_CLASS, BOTTOM_CLASS);  %BOTTOM_CLASS is changed automatically
TOP_CLASS = add_CHILD_snow(TOP_CLASS, class_list, stratigraphy_list);

CURRENT = TOP.NEXT;
while ~isequal(CURRENT, BOTTOM)
    plot(CURRENT.STATVAR.T, CURRENT.STATVAR.upperPos - cumsum(CURRENT.STATVAR.layerThick))
    hold on
    fprintf('Upper Position: %3.1f \nTemperature: %2.1f\n', CURRENT.STATVAR.upperPos, CURRENT.STATVAR.T(1))
    fprintf('Lower Position: %3.1f \nTemperature: %2.1f\n', CURRENT.STATVAR.lowerPos, CURRENT.STATVAR.T(end))
    CURRENT = CURRENT.NEXT;
end

%------ time integration ------------------
t = forcing.PARA.start_time;

%t is in days, timestep should also be in days
while t <= forcing.PARA.end_time

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
    timestep = 10; %3600; in days!
    while ~isequal(CURRENT, BOTTOM)
         timestep = min(timestep, get_timestep(CURRENT));
        CURRENT = CURRENT.NEXT;
    end
    %timestep = min(timestep, (out.OUTPUT_TIME-t).*day_sec);
    timestep = min(timestep, (out.OUTPUT_TIME-t)); %in days!!
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
    CURRENT = BOTTOM.PREVIOUS;
    while ~isequal(CURRENT, TOP)
        CURRENT = compute_diagnostic(CURRENT, forcing);
        CURRENT = CURRENT.PREVIOUS;
    end
    
    TOP_CLASS = TOP.NEXT; %TOP_CLASS and BOTOOM_CLASS for convenient access
    BOTTOM_CLASS = BOTTOM.PREVIOUS;


    %store the output according to the defined OUT clas
    out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number);
    
    %calculate new time
    t = t + timestep; %./day_sec;
                      
    if out.BREAK == 1
        break
    end
    
end


display_ground_SubseaPFwSalt(strcat('results', filesep, run_number, filesep, 'Results_', run_number))
ok = 1;




