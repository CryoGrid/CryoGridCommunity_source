classdef PROVIDER_YAML3D < PROVIDER_YAML
    
    
    methods
        

        function provider = assign_paths_yaml3d(provider, run_name, result_path, constant_file)
            provider = assign_paths_yaml(provider, run_name, result_path, constant_file);
        end
       
        
        function provider = read_const_yaml3d(provider)
            provider = read_const_yaml(provider);
        end
        
        function provider = read_parameters_yaml3d(provider)
            provider = read_parameters_yaml(provider);
        end
           
        %parallel runs
        function provider = update_parameter_file_yaml3d(provider, param_file_number)
            %change parameter file
            provider = update_parameter_file_yaml(provider, param_file_number);
        end
        
        function provider = update_run_name_yaml3d(provider, worker_number)
            %change run name and make new result directories if necessary
            provider = update_run_name_yaml(provider, worker_number);
        end
        
        
        
    end
end



