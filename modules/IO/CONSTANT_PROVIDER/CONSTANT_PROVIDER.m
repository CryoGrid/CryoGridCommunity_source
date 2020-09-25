function cprovider = CONSTANT_PROVIDER(config_path, init_format, const_file)

if isequal(init_format, 'EXCEL')
    const_file = [const_file '.xlsx'];
    cprovider = CONSTANT_PROVIDER_EXCEL(config_path, const_file);
    
elseif isequal(init_format, 'YAML')
    const_file = [const_file '.yml'];
    cprovider = CONSTANT_PROVIDER_YAML(config_path, const_file);
else
    cprovider = 0;
end