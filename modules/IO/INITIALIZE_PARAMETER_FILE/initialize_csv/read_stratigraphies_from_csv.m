function [ stratigraphy_list ] = read_stratigraphies_from_csv( path_to_configuration_file, available_classes_list )
% This function creates the different STRAT objects (namely STRAT_linear, STRAT_layers and STRAT_classes), reads the values from the corresponding csv
% configuration files and initializes the main variables. Also calls
% several functions to check the validity of the inputs.  

STRAT_classes_names = available_classes_list(contains(available_classes_list,'STRAT') == 1);

stratigraphy_list ={};
for i = 1:length(STRAT_classes_names)
    STRAT_class_name = STRAT_classes_names{i,1};
    t = readtable([path_to_configuration_file '\' STRAT_class_name '.csv']);
    strat_csv_inputs_check(t, STRAT_class_name, available_classes_list)
    class_handle = str2func(STRAT_class_name);
    strat = class_handle();
    strat = initialize_from_table(strat, t);
    stratigraphy_list = [stratigraphy_list; {strat}];
    stratigraphy_list{i,2} = 1;
end


end

