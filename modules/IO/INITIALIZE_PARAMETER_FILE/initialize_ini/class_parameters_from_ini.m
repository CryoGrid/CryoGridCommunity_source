
function class_list = class_parameters_from_ini( path_to_config_directory, available_class_names )
% This function creates the different ground and snow objects from the corresponding classes, reads the values from the ini
% configuration file and initializes the main variables. Also calls
% several functions to check the validity of the inputs.    

    %WARNING !!! Ini parser does not handle sections with the same name,
    %problematic to configure similar classes with different variables'
    %values. For now, each class only appears once in the ini file without
    %mentionning the index ... xml configuration files might be the solution in the future
	
    % Ini parser importing the parameters' values from the configuration file and storing them into a
    % structure
    ini = ini2struct([path_to_config_directory '\CLASS.ini']);
    % Checking the validity and range of the inputs from the ini file
    classes_ini_inputs_check(ini);
    
    % Extracting class names from the section headers of the configuration file
    ini_classes_names = fieldnames(ini);
    
    % Checking whether  the classes exist within the model's folder structure and creating the corresponding object from class cited in the
    % configuration file
	class_list ={};
	for i = 1:length(ini_classes_names)
        class_exists(ini_classes_names{i,1}, available_class_names);
		class_handle = str2func(ini_classes_names{i,1});
		class = class_handle();
		class = provide_variables(class);
    
		% Initialization with the configuration file parameters
		class.PARA = initialize_from_ini(class, class.PARA, ini_classes_names{i,1}, ini); 
		
		class_list = [class_list; {class}];
        class_list{i,2} = 1;
	end
       
end



