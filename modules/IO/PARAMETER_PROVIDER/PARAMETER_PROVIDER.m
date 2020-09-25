function pprovider = PARAMETER_PROVIDER(config_path, init_format, run_number)

if isequal(init_format, 'EXCEL')
    parameter_file = [run_number '.xlsx'];
    pprovider = PARAMETER_PROVIDER_EXCEL(config_path, parameter_file);
    
elseif isequal(init_format, 'YAML')
    parameter_file = [run_number '.yml'];
    pprovider = PARAMETER_PROVIDER_YAML(config_path, parameter_file);
else 
    pprovider = 0;
end