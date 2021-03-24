classdef PROVIDER < PROVIDER_EXCEL & PROVIDER_EXCEL3D & PROVIDER_MAT
    
    methods
        
        %possibly work with varargin to allow for different initialization paramters for other PROVIDER classes
        function provider = assign_paths(provider, init_format, run_name, result_path, constant_file)
            provider.PARA.init_format = init_format;
            if strcmp(init_format, 'EXCEL')
                provider = assign_paths_excel(provider, run_name, result_path, constant_file);
            elseif strcmp(init_format, 'EXCEL3D')
                provider = assign_paths_excel3d(provider, run_name, result_path, constant_file);
            elseif strcmp(init_format, 'MAT')
                provider = assign_paths_mat(provider, run_name, result_path);
            end
        end
        
        function provider = read_const(provider)
            if strcmp(provider.PARA.init_format, 'EXCEL')
                provider = read_const_excel(provider);
            elseif strcmp(provider.PARA.init_format, 'EXCEL3D')
                provider = read_const_excel3d(provider);
            end
        end

        
        function provider = read_parameters(provider)
            if strcmp(provider.PARA.init_format, 'EXCEL')
                provider = read_parameters_excel(provider);
            elseif strcmp(provider.PARA.init_format, 'EXCEL3D')
                provider = read_parameters_excel3d(provider);
            end
        end
        
       function provider = update_parameter_file(provider, param_file_number)
           if strcmp(provider.PARA.init_format, 'EXCEL3D')
               provider = update_parameter_file_excel3d(provider, param_file_number);
           end
       end
       
       function provider = update_run_name(provider, worker_number)
           if strcmp(provider.PARA.init_format, 'EXCEL3D')
               provider = update_run_name_excel3d(provider, worker_number);
           end
       end
        
    end
end

