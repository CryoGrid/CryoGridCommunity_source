function class_exists( config_class_name, available_classes_list )
% This function verifies whether the name of the classes within the ini
% configuration corresponds to an existing class of the model's folder
% structure.

if any(strcmp(available_classes_list, config_class_name)) == 0
    error('This class is not available. Modify the class name in the section header of the .ini file')
end
end

