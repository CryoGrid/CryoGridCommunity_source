function out_ini_inputs_check( ini_struct )
% This function verifies the validity of the potential variables' inputs
% extracted from the resulting structure of the ini parser.

ini_out_name = fieldnames(ini_struct);
ini_out_variables = fields(ini_struct.(ini_out_name{1,1}));

%Check the data type, extension, range, etc. of the different inputs from
%potentially existing variables

if ~isempty(find(strcmp(ini_out_variables, 'output_timestep') == 1))
    i1 = find(strcmp(ini_out_variables, 'output_timestep') == 1);
    validateattributes(ini_struct.(ini_out_name{1}).(ini_out_variables{i1,1}), {'numeric'}, {'positive'});
end

if ~isempty(find(strcmp(ini_out_variables, 'save_date') == 1))
    i2 = find(strcmp(ini_out_variables, 'save_date') == 1);
    is_valid_date_ddmm(ini_struct.(ini_out_name{1}).(ini_out_variables{i2,1}));
end

if ~isempty(find(strcmp(ini_out_variables, 'save_interval') == 1))
    i3 = find(strcmp(ini_out_variables, 'save_interval') == 1);
    validateattributes(ini_struct.(ini_out_name{1}).(ini_out_variables{i3,1}), {'numeric'}, {'positive'});
end

end