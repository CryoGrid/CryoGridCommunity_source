%% Initialize CryoGrid - required before accessing output data
% R. Zweigel June 2019
% Equivalent to first 40 lines of main.m
% Must be run while in the CryoGrid folder where this script is found

clear all
addpath(genpath('modules'))

run_number = 'test_excel'; % File where classes are defined
result_path = 'results/';
parameter_file_type = 'xlsx';
const_file = 'CONSTANTS_excel.xlsx';

%read information for all classes from Excel file and store it in cell
%arrays
parameter_info = read_excel2cell(['template.xlsx']);  %read excel file to cell array

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

clear all
