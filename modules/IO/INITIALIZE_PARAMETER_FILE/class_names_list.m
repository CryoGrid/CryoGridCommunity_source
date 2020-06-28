function [ available_classes ] = class_names_list( path_to_modules )
%This function lists all the existing and available classes within the working environment directory that can be
%used by the model to run.

% Lists all the modules and their features available within the folder
% 'modules'
module_names = dir(path_to_modules);
module_names(ismember({module_names.name},{'.','..'})) = [];

% Lists all files and their features, from each module folder
sub_module.files =[];
for i = 1:length(module_names)
    sub_module(i).files = dir([path_to_modules '\' module_names(i).name '\**\*.m']);     
end

% Lists and stores the names of all existing classes within the 'modules' folder 
available_class_names = [];
for i = 1:length(sub_module)
    for j = 1:length(sub_module(i).files)
        [pathstr,name,ext] = fileparts(sub_module(i).files(j).name);
        sub_module(i).files(j).file_name = name;
        sub_module(i).files(j).isclass = isobject(str2num(sub_module(i).files(j).file_name)) == 1;
        if sub_module(i).files(j).isclass == 1
            module_names(i).model_classes(j).class_names = sub_module(i).files(j).file_name;
        end
    end
    available_class_names = [available_class_names module_names(i).model_classes];
end

% Only stores class names while removing empty fields
available_classes =[];
for i = 1:length(available_class_names)
    if ischar(available_class_names(i).class_names) == 1
    available_classes =[available_classes; {available_class_names(i).class_names}];
    end
end

end

