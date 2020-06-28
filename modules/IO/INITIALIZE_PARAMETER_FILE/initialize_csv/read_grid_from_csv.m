function grid = read_grid_from_csv(path_to_configuration_file, available_classes_list)
% This function creates the grid object, reads the values from the csv
% configuration file and initializes the main variables. Also calls
% several functions to check the validity of the inputs.  

    % Extract GRID class name from available class list
    grid_class_name = available_classes_list(contains(available_classes_list,'GRID') == 1)
    grid_class_name = cell2mat(grid_class_name);
    
    % Creating the corresponding object from GRID class 
    class_handle = str2func(grid_class_name);
    grid = class_handle();
    
    % Converting csv configuration file in table and checking the validity
    % of the inputs
    t = readtable([path_to_configuration_file '\' grid_class_name '.csv']);
    grid_csv_inputs_check(t);
    % check the option uiimport to optimize code ? and save code as new
    % function ? or check csvread, but does not seem to handle strings as
    % headers
    grid = grid.initialize_from_table(t);
end
