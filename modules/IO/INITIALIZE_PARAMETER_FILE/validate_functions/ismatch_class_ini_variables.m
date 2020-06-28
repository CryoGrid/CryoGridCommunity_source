function ismatch_class_ini_variables( ini_variables, class_variables )
% This function checks whether the variable names from the configuration
% file, match the variable names of the class being initialized.

if (all(strcmp(ini_variables, class_variables)) == 0) | (length(ini_variables) ~= length(class_variables))
    check = horzcat(ini_variables, class_variables)
    error('The variable names from the class differ from the variable names from the configuration file. Please make sure the variable names of the configuration file match.')
end

end

