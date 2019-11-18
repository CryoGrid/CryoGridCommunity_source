% CryoGrid main file to be executed
% S.Westermann, Jan 2019

clear all
close all

mod = genpath('modules');
addpath(mod);


%results are stored in results/"run_number"
%this folder must contain an Excel spreadsheet called "run_number".xlsx deifining the properties of all classes requried for the run and the
%spreadsheet "CONSTANTS_excel.xlsx" containing all constants

run_number = 'test_excel';
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
% TOP_CLASS = add_CHILD_snow(TOP_CLASS, class_list, stratigraphy_list); % commented out for now
 
%------ time integration ------------------
day_sec = 24.*3600;

t = forcing.PARA.start_time; 

tic()
%t is in days, timestep should also be in days
while t <= forcing.PARA.end_time

    [forcing] = interpolate_forcing(t, forcing);
    
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
    CURRENT = get_boundary_condition_l(CURRENT, forcing);  %At this point, CURRENT is equal to BOTTOM_CLASS 
    %--------------------------


    
    %calculate spatial derivatives for every cell in the stratigraphy
    CURRENT = TOP.NEXT;
    while ~isequal(CURRENT, BOTTOM)
        CURRENT = get_derivatives_prognostic(CURRENT);
        CURRENT = CURRENT.NEXT;
    end
    
    %calculate minimum timestep required for all cells in days
    CURRENT = TOP.NEXT;
    timestep = 3600;
    while ~isequal(CURRENT, BOTTOM)
         timestep = min(timestep, get_timestep(CURRENT));
        CURRENT = CURRENT.NEXT;
    end
    timestep = min(timestep, (out.OUTPUT_TIME-t).*day_sec);
%     make sure to hit the output times!
     
        
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
    
    %calculate new time
    t = t + timestep./day_sec;
   
    %store the output according to the defined OUT class
    out = store_OUT(out, t, TOP_CLASS, BOTTOM, forcing, run_number);
    
end
toc()


% maxT = max(out.STRATIGRAPHY(:,2));

% fname = 'C:\Users\sstuenzi\Seafile\sync\AWI\programing\CryoGrid\coupling\20191021_CryoVeg_Git\Output';
% mkdir 'C:\Users\sstuenzi\Seafile\sync\AWI\programing\CryoGrid\coupling\20191021_CryoVeg_Git\Output\20191029'
% fname = 'C:\Users\sstuenzi\Seafile\sync\AWI\programing\CryoGrid\coupling\20191021_CryoVeg_Git\Output\20191029';
% saveas(figure(20), 'Ground_temp.fig');
% saveas(figure(1), fullfile(fname, 'Vegetation_temp'), 'fig');
% saveas(figure(2), fullfile(fname, 'Storage_heat_flux'), 'fig');
% saveas(figure(3), fullfile(fname, 'Sensible_heat_flux'), 'fig');
% saveas(figure(9), fullfile(fname, 'Photosnthesis'), 'fig');
% saveas(figure(4), fullfile(fname, 'Latent_heat_flux'), 'fig');
% saveas(figure(5), fullfile(fname, 'Net_radiation'), 'fig');
% saveas(figure(6), fullfile(fname, 'Wind_profile'), 'fig');
% saveas(figure(11), fullfile(fname, 'Sunny_shaded_leaves'), 'fig');


