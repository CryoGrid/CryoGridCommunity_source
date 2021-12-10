classdef PROVIDER_EXCEL3D < PROVIDER_EXCEL
    
    
    methods
        

        function provider = assign_paths_excel3d(provider, run_name, result_path, constant_file)
            provider = assign_paths_excel(provider, run_name, result_path, constant_file);
        end
       
        
        function provider = read_const_excel3d(provider)
            provider = read_const_excel(provider);
        end
        
        function provider = read_parameters_excel3d(provider)
            provider = read_parameters_excel(provider);
        end
           
        %parallel runs
        function provider = update_parameter_file_excel3d(provider, param_file_number)
            %change parameter file
            parameter_file = [provider.PARA.result_path provider.PARA.run_name '/' provider.PARA.run_name '_' num2str(param_file_number) '.xlsx'];
            provider.PARA.parameter_file = parameter_file;
        end
        
        function provider = update_run_name_excel3d(provider, worker_number)
            %change run name and make new result directories if necessary
            provider.PARA.run_name = [provider.PARA.run_name '_' num2str(worker_number)];
            if ~(exist([provider.PARA.result_path provider.PARA.run_name ])==7)
                mkdir([provider.PARA.result_path provider.PARA.run_name]);
            end
        end

        
    end
end

